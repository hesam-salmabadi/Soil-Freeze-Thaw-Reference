"""Site-year usability assessment (fit-readiness gate).

Decides, per freezing cycle, whether there is enough data to fit an SFCC. The
real criterion is **temperature-axis coverage**: a freezing cycle is usable when
its points span the transition band without a hole larger than ``max_temp_gap``
and it reaches at least ``min_reach``. Time-axis coverage (critical-window
coverage %, longest outage) is computed too, but only as diagnostic flags — a
gap only hurts if it erases part of the temperature transition.

Outputs:
    - the labeled frame with a boolean ``usable`` column added (True for rows in a
      usable freezing cycle; False for thawing / no-cycle / unusable rows),
    - a per-cycle report DataFrame (one row per freezing cycle), writable as a
      sidecar CSV via :func:`write_report`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..preprocess.cycles import FREEZING


def _assess_freezing_cycle(
    temps: np.ndarray,
    *,
    warm_edge: float,
    cold_edge: float,
    min_reach: float,
    max_temp_gap: float,
) -> dict:
    """Assess one freezing cycle's temperature-axis coverage."""
    t = np.sort(np.asarray(temps, dtype=float))
    t = t[~np.isnan(t)]
    if t.size == 0:
        return {
            "n_points": 0, "t_min": np.nan, "t_max": np.nan, "shallow": True,
            "max_temp_hole": np.inf, "usable": False, "reason": "no valid points",
        }

    t_min, t_max = float(t[0]), float(t[-1])
    reached = t_min <= min_reach

    # Check for holes from the warm edge down to the deepest temperature we care
    # about (the cycle's min, but no deeper than cold_edge).
    lo = max(cold_edge, t_min)
    inband = t[(t >= lo) & (t <= warm_edge)]
    checkpoints = np.unique(np.concatenate(([lo, warm_edge], inband)))
    gaps = np.diff(checkpoints)
    max_hole = float(gaps.max()) if gaps.size else 0.0

    coverage_ok = max_hole <= max_temp_gap
    usable = bool(coverage_ok and reached)

    reasons = []
    if not coverage_ok:
        reasons.append(f"temp hole {max_hole:.2f}C > {max_temp_gap}C")
    if not reached:
        reasons.append(f"shallow (min {t_min:.2f}C > {min_reach}C)")
    reason = "ok" if usable else "; ".join(reasons)

    return {
        "n_points": int(t.size), "t_min": t_min, "t_max": t_max,
        "shallow": not reached, "max_temp_hole": max_hole,
        "usable": usable, "reason": reason,
    }


def _time_diagnostics(
    year_index: pd.DatetimeIndex, year: int, cw_start: tuple[int, int], cw_end: tuple[int, int]
) -> dict:
    """Critical-window time coverage and longest outage (diagnostic only)."""
    start = pd.Timestamp(year=year, month=cw_start[0], day=cw_start[1])
    end = pd.Timestamp(year=year + 1, month=cw_end[0], day=cw_end[1])
    win = year_index[(year_index >= start) & (year_index < end)].sort_values()
    total_hours = (end - start) / pd.Timedelta("1h")

    if len(win) == 0:
        return {"time_coverage": 0.0, "longest_gap_days": (end - start) / pd.Timedelta("1D")}

    present_hours = win.floor("1h").nunique()
    coverage = min(1.0, present_hours / total_hours)

    edges = pd.DatetimeIndex([start, *win, end])
    longest_gap = edges.to_series().diff().max()
    return {
        "time_coverage": round(float(coverage), 4),
        "longest_gap_days": round(longest_gap / pd.Timedelta("1D"), 3),
    }


def assess_usability(
    labeled: pd.DataFrame,
    *,
    temp_col: str = "soil_temperature",
    value_col: str | None = None,
    warm_edge: float = 2.0,
    cold_edge: float = -2.0,
    min_reach: float = -1.0,
    max_temp_gap: float = 0.5,
    critical_window: tuple[tuple[int, int], tuple[int, int]] = ((10, 1), (3, 1)),
    site_id: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Add a ``usable`` column and return (labeled_with_usable, per_cycle_report).

    A freezing cycle is ``usable`` when its temperature points cover the transition
    band with no hole larger than ``max_temp_gap`` and it reaches ``min_reach``. If
    ``value_col`` is given, only rows with that column non-NaN count as data points.
    """
    out = labeled.copy()
    usable_col = np.zeros(len(out), dtype=bool)
    rows: list[dict] = []

    is_freezing = out["cycle_phase"] == FREEZING
    cw_start, cw_end = critical_window

    for (year, cid), grp in out[is_freezing].groupby(["freezing_year", "cycle_id"]):
        if year < 0 or cid < 0:
            continue
        valid = grp[temp_col].notna()
        if value_col is not None and value_col in grp.columns:
            valid &= grp[value_col].notna()
        temps = grp.loc[valid, temp_col].to_numpy()

        assessment = _assess_freezing_cycle(
            temps, warm_edge=warm_edge, cold_edge=cold_edge,
            min_reach=min_reach, max_temp_gap=max_temp_gap,
        )
        year_index = out.index[out["freezing_year"] == year]
        assessment.update(_time_diagnostics(year_index, int(year), cw_start, cw_end))

        if assessment["usable"]:
            usable_col[(out["freezing_year"] == year) & (out["cycle_id"] == cid)] = True

        row = {"freezing_year": int(year), "cycle_id": int(cid)}
        if site_id is not None:
            row = {"site_id": site_id, **row}
        row.update(assessment)
        rows.append(row)

    out["usable"] = usable_col
    columns = (
        (["site_id"] if site_id is not None else [])
        + ["freezing_year", "cycle_id", "n_points", "t_min", "t_max",
           "shallow", "max_temp_hole", "usable",
           "time_coverage", "longest_gap_days", "reason"]
    )
    report = pd.DataFrame(rows)
    if not report.empty:
        report = report[[c for c in columns if c in report.columns]]
    return out, report


def assess_with_config(labeled, config, *, temp_col, site_id=None):
    """Convenience wrapper that unpacks a :class:`softer.config.Config`."""
    u = config.usability
    cw = config.coverage.critical_window
    return assess_usability(
        labeled,
        temp_col=temp_col,
        value_col=u.value_col,
        warm_edge=u.warm_edge,
        cold_edge=u.cold_edge,
        min_reach=u.min_reach,
        max_temp_gap=u.max_temp_gap,
        critical_window=(cw.start, cw.end),
        site_id=site_id,
    )


def write_report(report: pd.DataFrame, path: str) -> None:
    """Write the per-cycle usability report to a sidecar CSV."""
    report.to_csv(path, index=False)
