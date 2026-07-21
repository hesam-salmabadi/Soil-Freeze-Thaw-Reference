"""Freeze/thaw cycle detection and labeling (Sect. 2.2).

Takes a QA'd per-site time series and adds per-row labels marking which
freeze/thaw cycle each observation belongs to. Only the timestamp and
soil-temperature columns drive the logic; every other column (raw sensor
output, effective permittivity, ...) rides along untouched. No interpolation is
performed — detection works on the available samples, so the method is
insensitive to occasional missing data.

Core method (unchanged from the author's original workflow):
    Each reading is classified against the deadband h (= 2*sigma_t):
        warm  (T > +h),  cold  (T < -h),  neutral (|T| <= h).
    A *freezing* cycle is a warm block followed by a cold block (the soil truly
    dropped below -h); a *thawing* cycle is a cold block followed by a warm block.
    Neutral near-zero readings are absorbed into the surrounding state, and legs
    shorter than ``min_duration`` are dissolved into their neighbour, so diurnal
    noise near 0 degC never spawns spurious cycles.

Per-year handling (Aug 1 -> next Aug 1 windows):
    - Normal winter (soil reaches warm before its first cold): full multi-cycle
      detection, keeping every freezing cycle.
    - Cold start (window begins already frozen, e.g. Aug at -4 degC): flagged
      ``start_frozen`` and split at the coldest point (start -> min = freezing,
      min -> end = thawing).
    - Never-freezing (never below -h) -> year_status = 1; always-frozen (never
      above +h) -> year_status = 2. Both use the coldest-point split.

Added columns:
    - ``freezing_year`` : Aug-start year of the window a row falls in (``-1`` if none).
    - ``cycle_id``      : 0-based per freezing year; ``-1`` outside any cycle.
    - ``cycle_phase``   : ``1`` = freezing, ``0`` = thawing, ``-1`` = no cycle.
    - ``year_status``   : ``0`` normal, ``1`` never-freezing, ``2`` always-frozen.

Per-year diagnostics (method, start_frozen, partial_year, n_cycles, ...) are
returned in ``df.attrs["years"]``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..sensors import get_sigma_t

# cycle_phase / cycle_id sentinels.
FREEZING = 1
THAWING = 0
NO_CYCLE = -1

# year_status values.
YEAR_NORMAL = 0
YEAR_NEVER_FROZEN = 1
YEAR_ALWAYS_FROZEN = 2

# Amplitude (degC) below which a window is treated as degenerate/flat.
_DEGENERATE_EPS = 1e-6


def _resolve_deadband(
    deadband_c: float | None, sigma_t: float | None, sensor: str | None
) -> float:
    """Resolve the deadband half-width h (degC).

    Priority: ``deadband_c`` > ``2 * sigma_t`` > ``2 * sigma_t(sensor)``.
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
    df: pd.DataFrame, temp_col: str, time_col: str | None
) -> pd.DataFrame:
    """Return a copy indexed by a sorted, de-duplicated DatetimeIndex (no interpolation)."""
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
    return work[~work.index.duplicated(keep="first")]


def _freezing_year(index: pd.DatetimeIndex, year_start: tuple[int, int]) -> np.ndarray:
    """Map each timestamp to the Aug-start year of its ``[start, next start)`` window."""
    month, day = year_start
    start_this = pd.to_datetime(
        {"year": index.year, "month": month, "day": day}
    ).to_numpy()
    before = index.to_numpy() < start_this
    years = index.year.to_numpy().copy()
    years[before] -= 1
    return years


def _split_at_gaps(idx: pd.DatetimeIndex, max_gap: pd.Timedelta | None) -> list[slice]:
    """Positional slices of ``idx`` split wherever the sample gap exceeds ``max_gap``."""
    n = len(idx)
    if n == 0:
        return []
    if max_gap is None or n == 1:
        return [slice(0, n)]
    deltas = np.diff(idx.to_numpy()).astype("timedelta64[ns]")
    breaks = np.nonzero(deltas > np.timedelta64(max_gap))[0] + 1
    bounds = [0, *breaks.tolist(), n]
    return [slice(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]


def _thermal_sign(valid_temp: pd.Series, h: float) -> pd.Series:
    """Warm (+1) / cold (-1) state per present reading; neutrals absorbed into neighbours."""
    raw = np.where(valid_temp > h, 1.0, np.where(valid_temp < -h, -1.0, np.nan))
    return pd.Series(raw, index=valid_temp.index).ffill().bfill()


def _block_anchors(valid_temp: pd.Series, sign: pd.Series) -> list[tuple[pd.Timestamp, int]]:
    """One extremum anchor per maximal same-sign block: peak if warm, trough if cold."""
    block = (sign != sign.shift()).cumsum()
    anchors: list[tuple[pd.Timestamp, int]] = []
    for _, idx in valid_temp.groupby(block).groups.items():
        cls = int(sign.loc[idx[0]])
        seg = valid_temp.loc[idx]
        anchors.append((seg.idxmax() if cls > 0 else seg.idxmin(), cls))
    return anchors


def _legs_from_anchors(
    anchors: list[tuple[pd.Timestamp, int]],
    seg_start: pd.Timestamp,
    seg_end: pd.Timestamp,
) -> list[list]:
    """Build [start, end, phase] legs between alternating anchors, extending shoulders."""
    legs: list[list] = []
    for i in range(len(anchors) - 1):
        t0, cls0 = anchors[i]
        t1, _ = anchors[i + 1]
        phase = FREEZING if cls0 > 0 else THAWING  # peak->trough cools = freezing
        legs.append([t0, t1, phase])
    if legs:
        legs[0][0] = seg_start  # extend first leg back to the segment start
        legs[-1][1] = seg_end   # extend last leg forward to the segment end
    return legs


def _minsplit_legs(
    valid_temp: pd.Series, seg_start: pd.Timestamp, seg_end: pd.Timestamp
) -> list[list]:
    """Coldest-point split: start -> min = freezing, min -> end = thawing."""
    if valid_temp.empty or (valid_temp.max() - valid_temp.min()) < _DEGENERATE_EPS:
        return []
    tmin = valid_temp.idxmin()
    legs: list[list] = []
    if tmin > seg_start:
        legs.append([seg_start, tmin, FREEZING])
    legs.append([tmin, seg_end, THAWING])
    return legs


def _absorb_short_legs(legs: list[list], min_duration: pd.Timedelta | None) -> list[list]:
    """Dissolve legs shorter than min_duration into their neighbour (merge, never delete-to-none)."""
    if min_duration is None or len(legs) <= 1:
        return legs
    legs = [list(leg) for leg in legs]
    while len(legs) > 1:
        durations = [t1 - t0 for t0, t1, _ in legs]
        short = [i for i, d in enumerate(durations) if d < min_duration]
        if not short:
            break
        i = min(short, key=lambda k: durations[k])
        if i == 0:
            legs[1][0] = legs[0][0]
            legs.pop(0)
        elif i == len(legs) - 1:
            legs[-2][1] = legs[-1][1]
            legs.pop()
        else:
            # Drop leg i; its two same-phase neighbours merge into one.
            legs[i - 1][1] = legs[i + 1][1]
            legs.pop(i + 1)
            legs.pop(i)
    return legs


def _segment_legs(
    seg_temp: pd.Series, method: str, h: float, min_duration: pd.Timedelta | None
) -> list[list]:
    """Return [start, end, phase] legs for one gap-free segment."""
    seg_temp = seg_temp.dropna()
    if len(seg_temp) < 2:
        return []
    seg_start, seg_end = seg_temp.index[0], seg_temp.index[-1]
    if method == "multicycle":
        sign = _thermal_sign(seg_temp, h)
        legs = _legs_from_anchors(_block_anchors(seg_temp, sign), seg_start, seg_end)
    else:  # "minsplit"
        legs = _minsplit_legs(seg_temp, seg_start, seg_end)
    return _absorb_short_legs(legs, min_duration)


def _classify_year(valid_temp: pd.Series, h: float) -> tuple[str, int, bool]:
    """Decide detection method and year_status for one freezing-year window.

    Returns (method, year_status, start_frozen).
    """
    below = valid_temp < -h
    above = valid_temp > h

    if not below.any():
        return "minsplit", YEAR_NEVER_FROZEN, False
    if not above.any():
        return "minsplit", YEAR_ALWAYS_FROZEN, False

    first_cold_time = valid_temp.index[below.to_numpy()][0]
    warm_before_cold = bool(above[valid_temp.index < first_cold_time].any())
    if warm_before_cold:
        return "multicycle", YEAR_NORMAL, False
    return "minsplit", YEAR_NORMAL, True  # cold start


def _assign_legs(index: pd.DatetimeIndex, legs: list[list]) -> tuple[np.ndarray, np.ndarray]:
    """Assign cycle_id / cycle_phase to every row by the leg its timestamp falls in."""
    cycle_id = np.full(len(index), NO_CYCLE, dtype=int)
    cycle_phase = np.full(len(index), NO_CYCLE, dtype=int)
    for i, (t0, t1, phase) in enumerate(legs):
        if i < len(legs) - 1:
            mask = (index >= t0) & (index < t1)
        else:
            mask = (index >= t0) & (index <= t1)
        cycle_id[mask] = i
        cycle_phase[mask] = phase
    return cycle_id, cycle_phase


def _is_partial_year(idx: pd.DatetimeIndex, year: int, year_start: tuple[int, int]) -> bool:
    """True if the window's data does not span most of the Aug->Aug year."""
    month, day = year_start
    start = pd.Timestamp(year=year, month=month, day=day)
    end = pd.Timestamp(year=year + 1, month=month, day=day)
    return (idx.max() - idx.min()) < (end - start) * 0.9


def label_freeze_thaw_cycles(
    df: pd.DataFrame,
    *,
    temp_col: str = "soil_temperature",
    time_col: str | None = None,
    sigma_t: float | None = None,
    deadband_c: float | None = None,
    sensor: str | None = None,
    year_start: tuple[int, int] = (8, 1),
    max_gap: str | pd.Timedelta | None = "2D",
    min_duration: str | pd.Timedelta | None = "2D",
) -> pd.DataFrame:
    """Label each row with its freeze/thaw cycle. See module docstring for the method.

    Parameters
    ----------
    df, temp_col, time_col:
        Input frame and the columns used. Non-temperature columns pass through.
    sigma_t, deadband_c, sensor:
        Set the deadband half-width h (priority: deadband_c > 2*sigma_t >
        2*sigma_t(sensor)). In production, sigma_t comes from per-site metadata.
    year_start:
        ``(month, day)`` start of the freezing year (default Aug 1).
    max_gap:
        A cycle never spans a data gap larger than this (``None`` disables).
    min_duration:
        Freeze/thaw legs shorter than this are absorbed into their neighbour.

    Returns
    -------
    pandas.DataFrame
        Sorted, datetime-indexed copy of ``df`` with ``freezing_year``, ``cycle_id``,
        ``cycle_phase`` and ``year_status`` columns added. ``df.attrs["years"]``
        holds per-year diagnostics.
    """
    h = _resolve_deadband(deadband_c, sigma_t, sensor)
    max_gap_td = pd.Timedelta(max_gap) if max_gap is not None else None
    min_dur_td = pd.Timedelta(min_duration) if min_duration is not None else None
    work = _as_sorted_datetime_frame(df, temp_col, time_col)

    n = len(work)
    freezing_year = np.full(n, NO_CYCLE, dtype=int)
    cycle_id = np.full(n, NO_CYCLE, dtype=int)
    cycle_phase = np.full(n, NO_CYCLE, dtype=int)
    year_status = np.full(n, NO_CYCLE, dtype=int)
    diagnostics: dict[int, dict] = {}

    if n == 0:
        for name, arr in (
            ("freezing_year", freezing_year),
            ("cycle_id", cycle_id),
            ("cycle_phase", cycle_phase),
            ("year_status", year_status),
        ):
            work[name] = arr
        work.attrs["years"] = diagnostics
        return work

    years = _freezing_year(work.index, year_start)

    for year in np.unique(years):
        ypos = np.nonzero(years == year)[0]
        y_index = work.index[ypos]
        y_temp = work[temp_col].iloc[ypos]
        valid = y_temp.dropna()

        if len(valid) < 2:
            freezing_year[ypos] = year
            diagnostics[int(year)] = {
                "method": "none",
                "year_status": None,
                "start_frozen": False,
                "partial_year": True,
                "n_cycles": 0,
            }
            continue

        method, status, start_frozen = _classify_year(valid, h)

        legs_all: list[list] = []
        for seg in _split_at_gaps(y_index, max_gap_td):
            legs_all.extend(_segment_legs(y_temp.iloc[seg], method, h, min_dur_td))

        y_cid, y_phase = _assign_legs(y_index, legs_all)
        freezing_year[ypos] = year
        cycle_id[ypos] = y_cid
        cycle_phase[ypos] = y_phase
        year_status[ypos] = status

        diagnostics[int(year)] = {
            "method": method,
            "year_status": status,
            "start_frozen": start_frozen,
            "never_frozen": status == YEAR_NEVER_FROZEN,
            "always_frozen": status == YEAR_ALWAYS_FROZEN,
            "partial_year": bool(_is_partial_year(y_index, int(year), year_start)),
            "n_cycles": len(legs_all),
        }

    work["freezing_year"] = freezing_year
    work["cycle_id"] = cycle_id
    work["cycle_phase"] = cycle_phase
    work["year_status"] = year_status
    work.attrs["years"] = diagnostics
    return work


def label_cycles_csv(in_path: str, out_path: str, **kwargs) -> pd.DataFrame:
    """Read a site CSV, label its cycles, write the result, and return the frame."""
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
