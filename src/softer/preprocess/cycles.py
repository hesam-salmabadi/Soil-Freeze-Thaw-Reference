"""Freeze/thaw cycle detection and labeling (Sect. 2.2).

Takes a QA'd per-site time series (produced by a separate upstream step) and adds
per-row labels marking which freeze/thaw cycle each observation belongs to. Only
the timestamp and soil-temperature columns drive the logic; every other column
(raw sensor output, effective permittivity, ...) rides along untouched.

Canonical logic (zero-crossing method):
    - A *freezing* cycle starts with soil temperature above 0 degC and declines
      continuously to a negative minimum.
    - A *thawing* cycle is the reverse: from that minimum up to above 0 degC.
    - Fluctuations within +/-2*sigma_t of 0 degC (sigma_t = instrument-specific
      temperature uncertainty) are ignored, so small wiggles near 0 do not spawn
      spurious cycles. This is the deadband.

This is a clean, row-preserving, side-effect-free reimplementation of the intent
of the original ``freeze_thaw_cycle_extractor_function_v2.process_cycles`` (which
dropped deadband rows and wrote one CSV per cycle). Plotting is kept out of the
compute path; use :func:`plot_cycles` for optional QA.

Added columns:
    - ``cycle_id``    : 0-based integer per detected cycle; ``-1`` outside any cycle.
    - ``cycle_phase`` : ``1`` = freezing, ``0`` = thawing, ``-1`` = no cycle.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..sensors import get_sigma_t

# Label constants.
FREEZING = 1
THAWING = 0
NO_CYCLE = -1


def _resolve_deadband(
    deadband_c: float | None,
    sigma_t: float | None,
    sensor: str | None,
) -> float:
    """Resolve the deadband half-width h (degC) from the provided inputs.

    Priority: explicit ``deadband_c`` > ``2 * sigma_t`` > ``2 * sigma_t(sensor)``.
    """
    if deadband_c is not None:
        h = float(deadband_c)
    elif sigma_t is not None:
        h = 2.0 * float(sigma_t)
    elif sensor is not None:
        h = 2.0 * get_sigma_t(sensor)
    else:
        raise ValueError(
            "Provide one of deadband_c, sigma_t, or sensor to set the "
            "near-zero deadband used for cycle detection."
        )
    if h < 0:
        raise ValueError(f"Deadband half-width must be non-negative, got {h}.")
    return h


def _as_sorted_datetime_frame(
    df: pd.DataFrame,
    temp_col: str,
    time_col: str | None,
    freq: str | None,
) -> pd.DataFrame:
    """Return a copy indexed by a sorted DatetimeIndex, optionally on a uniform grid."""
    work = df.copy()

    if time_col is not None:
        work.index = pd.to_datetime(work[time_col])
    elif not isinstance(work.index, pd.DatetimeIndex):
        raise ValueError(
            "Provide time_col, or pass a frame already indexed by a DatetimeIndex."
        )
    work = work.sort_index()

    if temp_col not in work.columns:
        raise KeyError(f"Temperature column {temp_col!r} not found in the frame.")

    if freq is not None:
        full = pd.date_range(work.index.min(), work.index.max(), freq=freq)
        work = work.reindex(full)
        num_cols = work.select_dtypes(include="number").columns
        obj_cols = work.select_dtypes(exclude="number").columns
        work[num_cols] = work[num_cols].interpolate(method="time")
        work[obj_cols] = work[obj_cols].ffill()

    return work


def _thermal_state(temp: pd.Series, h: float) -> pd.Series:
    """Classify each row as warm (+1) / cold (-1), absorbing near-zero neutrals.

    Rows within the deadband (|T| <= h) are neutral; they are filled with the
    surrounding non-neutral sign so a run of small near-zero wiggles cannot open
    a new cycle.
    """
    filled = temp.astype(float).interpolate(method="time").bfill().ffill()
    raw = np.where(filled > h, 1.0, np.where(filled < -h, -1.0, np.nan))
    sign = pd.Series(raw, index=temp.index).ffill().bfill()
    return sign


def _cycle_boundaries(temp: pd.Series, sign: pd.Series) -> list[tuple[int, int]]:
    """Return alternating extrema anchors as ``(positional_index, class)`` pairs.

    Each maximal same-sign block contributes one anchor: its temperature maximum
    if warm, its minimum if cold. Consecutive anchors therefore alternate
    peak/trough and delimit the cycles.
    """
    block = (sign != sign.shift()).cumsum()
    anchors: list[tuple[int, int]] = []
    for _, positions in block.groupby(block).groups.items():
        cls = int(sign.loc[positions[0]])
        segment = temp.loc[positions]
        anchor_label = segment.idxmax() if cls > 0 else segment.idxmin()
        anchors.append((temp.index.get_loc(anchor_label), cls))
    return anchors


def _assign_labels(n_rows: int, anchors: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    """Fill cycle_id / cycle_phase arrays from the ordered extrema anchors."""
    cycle_id = np.full(n_rows, NO_CYCLE, dtype=int)
    cycle_phase = np.full(n_rows, NO_CYCLE, dtype=int)

    for i in range(len(anchors) - 1):
        start, cls = anchors[i]
        end = anchors[i + 1][0]
        # The shared trough/peak row belongs to the following cycle, except for
        # the final cycle which includes its closing anchor.
        stop = end + 1 if i == len(anchors) - 2 else end
        phase = FREEZING if cls > 0 else THAWING
        cycle_id[start:stop] = i
        cycle_phase[start:stop] = phase

    return cycle_id, cycle_phase


def label_freeze_thaw_cycles(
    df: pd.DataFrame,
    *,
    temp_col: str = "soil_temperature",
    time_col: str | None = None,
    sigma_t: float | None = None,
    deadband_c: float | None = None,
    sensor: str | None = None,
    freq: str | None = None,
) -> pd.DataFrame:
    """Label each row with its freeze/thaw cycle.

    Parameters
    ----------
    df:
        Per-site time series. Must contain ``temp_col`` and either be indexed by a
        DatetimeIndex or carry ``time_col``.
    temp_col:
        Name of the soil-temperature column [degC].
    time_col:
        Name of the timestamp column. If ``None``, ``df`` must already be indexed
        by a DatetimeIndex.
    sigma_t, deadband_c, sensor:
        Set the near-zero deadband half-width ``h``. Priority: ``deadband_c`` >
        ``2 * sigma_t`` > ``2 * sigma_t`` looked up for ``sensor``.
    freq:
        Optional pandas offset (e.g. ``"1H"``). If given, the series is reindexed
        onto that uniform grid (numeric columns interpolated, object columns
        forward-filled) before detection. Off by default — the input is assumed
        already regular/QA'd.

    Returns
    -------
    pandas.DataFrame
        A sorted, datetime-indexed copy of ``df`` with integer ``cycle_id`` and
        ``cycle_phase`` columns added. ``df.attrs["never_frozen"]`` is set when the
        temperature never drops below ``-h`` (no freezing cycle exists).
    """
    h = _resolve_deadband(deadband_c, sigma_t, sensor)
    work = _as_sorted_datetime_frame(df, temp_col, time_col, freq)

    temp = work[temp_col].astype(float)
    temp_filled = temp.interpolate(method="time").bfill().ffill()

    never_frozen = not bool((temp_filled < -h).any())
    if never_frozen:
        work["cycle_id"] = NO_CYCLE
        work["cycle_phase"] = NO_CYCLE
        work.attrs["never_frozen"] = True
        return work

    sign = _thermal_state(temp, h)
    anchors = _cycle_boundaries(temp_filled, sign)
    cycle_id, cycle_phase = _assign_labels(len(work), anchors)

    work["cycle_id"] = cycle_id
    work["cycle_phase"] = cycle_phase
    work.attrs["never_frozen"] = False
    return work


def label_cycles_csv(in_path: str, out_path: str, **kwargs) -> pd.DataFrame:
    """Read a site CSV, label its cycles, write the result, and return the frame.

    ``kwargs`` are forwarded to :func:`label_freeze_thaw_cycles`. The timestamp
    column is inferred from ``time_col`` (if given) else assumed to be the first
    column / index of the CSV.
    """
    time_col = kwargs.get("time_col")
    if time_col is not None:
        frame = pd.read_csv(in_path)
    else:
        frame = pd.read_csv(in_path, index_col=0, parse_dates=True)
    labeled = label_freeze_thaw_cycles(frame, **kwargs)
    labeled.to_csv(out_path)
    return labeled


def plot_cycles(
    labeled: pd.DataFrame,
    *,
    temp_col: str = "soil_temperature",
    moisture_col: str | None = None,
):
    """Optional QA plot: temperature with freezing cycles shaded. Not used in the pipeline.

    Imports matplotlib lazily so the compute path has no plotting dependency.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(labeled.index, labeled[temp_col], color="gray", linewidth=1, label=temp_col)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=0.8)

    if "cycle_phase" in labeled.columns:
        freezing = labeled["cycle_phase"] == FREEZING
        # Shade contiguous freezing spans.
        span = (freezing != freezing.shift()).cumsum()
        for _, idx in labeled[freezing].groupby(span[freezing]).groups.items():
            ax.axvspan(idx[0], idx[-1], color="tab:blue", alpha=0.2)

    if moisture_col is not None and moisture_col in labeled.columns:
        ax2 = ax.twinx()
        ax2.plot(labeled.index, labeled[moisture_col], color="tab:purple", alpha=0.6)
        ax2.set_ylabel(moisture_col, color="tab:purple")

    ax.set_xlabel("Time")
    ax.set_ylabel(temp_col)
    ax.set_title("Freeze/thaw cycles (freezing shaded)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig
