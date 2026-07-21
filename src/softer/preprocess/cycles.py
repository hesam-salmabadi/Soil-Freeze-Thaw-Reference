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
      temperature uncertainty) are ignored (the deadband), so small wiggles near
      0 do not spawn spurious cycles.

Automation concerns handled here:
    - **Freezing years**: the record is split into ``[Aug 1, next Aug 1)`` windows
      and detection runs independently per window, so an end-of-winter thaw never
      merges with the next winter's freeze. ``cycle_id`` resets per year.
    - **Never-frozen / always-frozen fallback**: when a window has no complete
      warm->cold->warm cycle (temperature never drops below -h, or never rises
      above +h), a single cycle is still assigned as freezing = time(max T) ->
      time(min T), thawing = time(min T) -> window end.
    - **Data gaps**: cycles never span a time gap larger than ``max_gap``.
    - **Degenerate windows** (flat / too few points) are left unlabeled, flagged.

This is a clean, row-preserving, side-effect-free reimplementation of the intent
of the original ``freeze_thaw_cycle_extractor_function_v2.process_cycles``.
Plotting is kept out of the compute path; use :func:`plot_cycles` for optional QA.

Added columns:
    - ``freezing_year`` : Aug-start year of the window a row falls in (``-1`` if none).
    - ``cycle_id``      : 0-based per freezing year; ``-1`` outside any cycle.
    - ``cycle_phase``   : ``1`` = freezing, ``0`` = thawing, ``-1`` = no cycle.

Per-year diagnostics are returned in ``df.attrs["years"]``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..sensors import get_sigma_t

# Label constants.
FREEZING = 1
THAWING = 0
NO_CYCLE = -1

# Amplitude (degC) below which a window is considered degenerate/flat.
_DEGENERATE_EPS = 1e-6


def _resolve_deadband(
    deadband_c: float | None,
    sigma_t: float | None,
    sensor: str | None,
) -> float:
    """Resolve the deadband half-width h (degC).

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
    """Return a copy indexed by a sorted, de-duplicated DatetimeIndex."""
    work = df.copy()

    if time_col is not None:
        work.index = pd.to_datetime(work[time_col])
    elif not isinstance(work.index, pd.DatetimeIndex):
        raise ValueError(
            "Provide time_col, or pass a frame already indexed by a DatetimeIndex."
        )
    else:
        work.index = pd.to_datetime(work.index)

    if temp_col not in work.columns:
        raise KeyError(f"Temperature column {temp_col!r} not found in the frame.")

    work = work.sort_index()
    work = work[~work.index.duplicated(keep="first")]

    if freq is not None:
        full = pd.date_range(work.index.min(), work.index.max(), freq=freq)
        work = work.reindex(full)
        num_cols = work.select_dtypes(include="number").columns
        obj_cols = work.select_dtypes(exclude="number").columns
        work[num_cols] = work[num_cols].interpolate(method="time")
        work[obj_cols] = work[obj_cols].ffill()

    return work


def _freezing_year(index: pd.DatetimeIndex, year_start: tuple[int, int]) -> np.ndarray:
    """Map each timestamp to the Aug-start year of its ``[start, next start)`` window."""
    month, day = year_start
    start_this = pd.to_datetime(
        {"year": index.year, "month": month, "day": day}
    ).to_numpy()
    # Before the window start, a row belongs to the previous window's year.
    before = index.to_numpy() < start_this
    years = index.year.to_numpy().copy()
    years[before] -= 1
    return years


def _split_at_gaps(idx: pd.DatetimeIndex, max_gap: pd.Timedelta | None) -> list[slice]:
    """Return positional slices of ``idx`` split wherever the gap exceeds ``max_gap``."""
    n = len(idx)
    if n == 0:
        return []
    if max_gap is None or n == 1:
        return [slice(0, n)]
    deltas = np.diff(idx.to_numpy()).astype("timedelta64[ns]")
    breaks = np.nonzero(deltas > np.timedelta64(max_gap))[0] + 1
    bounds = [0, *breaks.tolist(), n]
    return [slice(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]


def _thermal_state(temp: pd.Series, h: float) -> pd.Series:
    """Classify each row as warm (+1) / cold (-1), absorbing near-zero neutrals."""
    filled = temp.astype(float).interpolate(method="time").bfill().ffill()
    raw = np.where(filled > h, 1.0, np.where(filled < -h, -1.0, np.nan))
    return pd.Series(raw, index=temp.index).ffill().bfill()


def _cycle_boundaries(temp: pd.Series, sign: pd.Series) -> list[tuple[int, int]]:
    """Return alternating extrema anchors as ``(positional_index, class)`` pairs.

    Positional indices are relative to ``temp`` (a single continuous segment).
    """
    if sign.isna().all():
        return []
    block = (sign != sign.shift()).cumsum()
    anchors: list[tuple[int, int]] = []
    for _, positions in block.groupby(block).groups.items():
        cls = int(sign.loc[positions[0]])
        segment = temp.loc[positions]
        anchor_label = segment.idxmax() if cls > 0 else segment.idxmin()
        anchors.append((temp.index.get_loc(anchor_label), cls))
    return anchors


def _fallback_boundaries(temp: pd.Series) -> list[tuple[int, int]]:
    """Max->min freezing then min->end thawing, for windows with no complete cycle.

    Returns anchors ``[(max_pos, warm), (min_pos, cold), (last_pos, warm)]`` so the
    generic assignment produces freezing (max->min) and thawing (min->end). Only
    the portion from the max onward is treated as a cycle; earlier rows stay -1.
    """
    filled = temp.astype(float).interpolate(method="time").bfill().ffill()
    if filled.empty or (filled.max() - filled.min()) < _DEGENERATE_EPS:
        return []
    max_pos = int(filled.to_numpy().argmax())
    min_pos = int(filled.to_numpy().argmin())
    if max_pos <= min_pos:
        # Normal freezing-year shape: warm start cooling to a later minimum.
        anchors = [(max_pos, 1), (min_pos, -1)]
        last = len(filled) - 1
        if last > min_pos:
            anchors.append((last, 1))  # trough -> end = thawing
        return anchors
    # Record starts cold and warms: treat min->max as a thawing leg.
    return [(min_pos, -1), (max_pos, 1)]


def _assign_labels(
    n_rows: int, anchors: list[tuple[int, int]]
) -> tuple[np.ndarray, np.ndarray]:
    """Fill cycle_id / cycle_phase arrays from ordered extrema anchors (segment-local)."""
    cycle_id = np.full(n_rows, NO_CYCLE, dtype=int)
    cycle_phase = np.full(n_rows, NO_CYCLE, dtype=int)

    for i in range(len(anchors) - 1):
        start, cls = anchors[i]
        end = anchors[i + 1][0]
        stop = end + 1 if i == len(anchors) - 2 else end
        cycle_id[start:stop] = i
        cycle_phase[start:stop] = FREEZING if cls > 0 else THAWING

    return cycle_id, cycle_phase


def _apply_min_filters(
    temp: pd.Series,
    cycle_id: np.ndarray,
    cycle_phase: np.ndarray,
    min_amplitude: float | None,
    min_duration: pd.Timedelta | None,
) -> None:
    """In place: drop cycles below min amplitude/duration back to NO_CYCLE."""
    if min_amplitude is None and min_duration is None:
        return
    for cid in np.unique(cycle_id):
        if cid < 0:
            continue
        mask = cycle_id == cid
        seg = temp.iloc[np.nonzero(mask)[0]]
        drop = False
        if min_amplitude is not None and (seg.max() - seg.min()) < min_amplitude:
            drop = True
        if min_duration is not None and (seg.index[-1] - seg.index[0]) < min_duration:
            drop = True
        if drop:
            cycle_id[mask] = NO_CYCLE
            cycle_phase[mask] = NO_CYCLE


def label_freeze_thaw_cycles(
    df: pd.DataFrame,
    *,
    temp_col: str = "soil_temperature",
    time_col: str | None = None,
    sigma_t: float | None = None,
    deadband_c: float | None = None,
    sensor: str | None = None,
    freq: str | None = None,
    year_start: tuple[int, int] = (8, 1),
    max_gap: str | pd.Timedelta | None = "2D",
    min_amplitude: float | None = None,
    min_duration: str | pd.Timedelta | None = None,
) -> pd.DataFrame:
    """Label each row with its freeze/thaw cycle.

    See module docstring for the algorithm and the added columns. Detection runs
    independently within each ``[Aug 1, next Aug 1)`` freezing-year window and
    within each gap-free segment of that window.

    Parameters
    ----------
    df, temp_col, time_col:
        Input frame and the columns used. Non-temperature columns pass through.
    sigma_t, deadband_c, sensor:
        Set the near-zero deadband half-width ``h`` (priority: deadband_c >
        2*sigma_t > 2*sigma_t(sensor)).
    freq:
        Optional pandas offset to reindex onto a uniform grid before detection.
    year_start:
        ``(month, day)`` start of the freezing year (default ``(8, 1)`` = Aug 1).
    max_gap:
        Cycles never span a time gap larger than this (``None`` disables).
    min_amplitude, min_duration:
        Optional filters dropping trivially small cycles (default off).

    Returns
    -------
    pandas.DataFrame
        Sorted, datetime-indexed copy of ``df`` with ``freezing_year``, ``cycle_id``
        and ``cycle_phase`` columns added. ``df.attrs["years"]`` holds per-year
        diagnostics (``partial_year``, ``detection``, ``never_frozen``, ``n_cycles``).
    """
    h = _resolve_deadband(deadband_c, sigma_t, sensor)
    max_gap_td = pd.Timedelta(max_gap) if max_gap is not None else None
    min_dur_td = pd.Timedelta(min_duration) if min_duration is not None else None
    work = _as_sorted_datetime_frame(df, temp_col, time_col, freq)

    n = len(work)
    freezing_year = np.full(n, NO_CYCLE, dtype=int)
    cycle_id = np.full(n, NO_CYCLE, dtype=int)
    cycle_phase = np.full(n, NO_CYCLE, dtype=int)
    diagnostics: dict[int, dict] = {}

    if n == 0:
        work["freezing_year"] = freezing_year
        work["cycle_id"] = cycle_id
        work["cycle_phase"] = cycle_phase
        work.attrs["years"] = diagnostics
        return work

    years = _freezing_year(work.index, year_start)

    for year in np.unique(years):
        ymask = years == year
        ypos = np.nonzero(ymask)[0]
        y_temp = work[temp_col].iloc[ypos]

        never_frozen = not bool((y_temp < -h).any())
        never_thawed = not bool((y_temp > h).any())
        next_id = 0
        used_fallback = False

        for seg in _split_at_gaps(y_temp.index, max_gap_td):
            seg_pos = ypos[seg]
            seg_temp = y_temp.iloc[seg]
            if len(seg_temp) < 2:
                continue

            sign = _thermal_state(seg_temp, h)
            anchors = _cycle_boundaries(seg_temp, sign)
            if len(anchors) < 2:
                anchors = _fallback_boundaries(seg_temp)
                if anchors:
                    used_fallback = True
            if len(anchors) < 2:
                continue

            seg_id, seg_phase = _assign_labels(len(seg_temp), anchors)
            _apply_min_filters(seg_temp, seg_id, seg_phase, min_amplitude, min_dur_td)

            valid = seg_id >= 0
            seg_id_shifted = np.where(valid, seg_id + next_id, NO_CYCLE)
            cycle_id[seg_pos] = seg_id_shifted
            cycle_phase[seg_pos] = seg_phase
            if valid.any():
                next_id = int(seg_id_shifted[valid].max()) + 1

        freezing_year[ypos] = year
        diagnostics[int(year)] = {
            "partial_year": bool(_is_partial_year(work.index[ypos], int(year), year_start)),
            "detection": "fallback" if used_fallback else "normal",
            "never_frozen": never_frozen,
            "never_thawed": never_thawed,
            "n_cycles": next_id,
        }

    work["freezing_year"] = freezing_year
    work["cycle_id"] = cycle_id
    work["cycle_phase"] = cycle_phase
    work.attrs["years"] = diagnostics
    return work


def _is_partial_year(
    idx: pd.DatetimeIndex, year: int, year_start: tuple[int, int]
) -> bool:
    """True if the window's data does not span most of the Aug->Aug year."""
    month, day = year_start
    start = pd.Timestamp(year=year, month=month, day=day)
    end = pd.Timestamp(year=year + 1, month=month, day=day)
    covered = idx.max() - idx.min()
    return covered < (end - start) * 0.9


def label_cycles_csv(in_path: str, out_path: str, **kwargs) -> pd.DataFrame:
    """Read a site CSV, label its cycles, write the result, and return the frame.

    ``kwargs`` are forwarded to :func:`label_freeze_thaw_cycles`.
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
    """Optional QA plot: temperature with freezing cycles shaded. Not used in the pipeline."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(labeled.index, labeled[temp_col], color="gray", linewidth=1, label=temp_col)
    ax.axhline(0.0, color="black", linestyle="--", linewidth=0.8)

    if "cycle_phase" in labeled.columns:
        freezing = labeled["cycle_phase"] == FREEZING
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
