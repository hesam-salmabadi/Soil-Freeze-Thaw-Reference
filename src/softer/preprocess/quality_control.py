"""Sequential QA/QC for near-surface soil temperature (design: ``documents/manuals/soil_temperature_qaqc_design.md``).

Automated, non-destructive quality control for multi-depth soil temperature time
series. Raw values are never modified; each test writes an ordinal flag to its own
companion column, rolled up only at the final aggregation stage so provenance
survives (a physically-valid zero-curtain flagged by persistence can be reviewed
rather than silently discarded).

Flag schema (QARTOD-ordinal)::

    1 PASS   2 NOT_EVALUATED   3 SUSPECT   4 FAIL   9 MISSING

``max`` over per-test flags gives the summary flag with the correct precedence
(9 > 4 > 3 > 2 > 1): a missing datum dominates, then FAIL, SUSPECT,
NOT_EVALUATED, PASS.

Processing order (each stage depth-aware; freeze-thaw logic in Sect. 6 is
cross-cutting and consumed by T4/T5)::

    T0  Structural           (missing, non-finite, gap segmentation)
    T1  Gross / sensor range (static physical bounds)
    T2a Climatology range    (per-depth per-DOY envelope, multi-year)
    T2b Ancillary consistency(ERA5 / land-surface-model residual; SUSPECT-capped)
    T3  Step / rate-of-change(per-hour derivative)
    T4  Spike                (Hampel/MAD, noise-floored)
    T5  Persistence          (T5a exact-repeat  +  T5b low-variance, FT-gated)
    T6  Profile consistency  (>=2 depths: amplitude / phase-lag / gradient)
    T7  Aggregation          (summary flag + quality metrics)

All window lengths and rate limits are specified in **physical time** (hours,
degC/hour) and converted to sample counts per site using the native interval, so
one config is valid across 30-min / hourly / 3-hourly sites. Shape/temporal tests
(T3-T6) run on contiguous segments and never span a gap.

Entry point: :func:`apply_qc`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Flag values (QARTOD-ordinal). Chosen so that ``max`` yields the right winner.
# --------------------------------------------------------------------------- #
PASS = 1
NOT_EVALUATED = 2
SUSPECT = 3
FAIL = 4
MISSING = 9

FLAG_VALUES = (1, 2, 3, 4, 9)
FLAG_MEANINGS = "pass not_evaluated suspect fail missing"

# Depth-band identifiers.
SURF = "SURF"  # ~0 cm skin/surface
SHAL = "SHAL"  # ~5 cm primary channel
SUB = "SUB"    # ~10 cm smoothest


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BandParams:
    """Per depth-band thresholds. Defaults are the design-doc appendix starting
    values (tune against each station's climatology and sensor noise floor)."""

    gross_min: float
    gross_max: float
    roc_suspect: float       # degC/hr
    roc_fail: float          # degC/hr
    spike_win_hours: float
    n_spike_suspect: float   # x robust sigma
    n_spike_fail: float      # x robust sigma
    sigma_noise: float       # degC; spike-scale noise floor (mandatory, Sect. 6)


_DEFAULT_BANDS: dict[str, BandParams] = {
    SURF: BandParams(-45.0, 70.0, 8.0, 12.0, 3.0, 4.0, 6.0, 0.10),
    SHAL: BandParams(-40.0, 55.0, 4.0, 6.0, 3.0, 3.0, 5.0, 0.08),
    SUB:  BandParams(-40.0, 50.0, 3.0, 5.0, 3.0, 3.0, 4.0, 0.05),
}

# Nominal depth (cm) -> band. Extend for other reported depths as needed.
_DEFAULT_DEPTH_TO_BAND: dict[str, str] = {"0cm": SURF, "5cm": SHAL, "10cm": SUB}


@dataclass(frozen=True)
class QCConfig:
    """All QC thresholds. Constant across sites; per-site facts (noise floor,
    delta_pc, present depths) should override before production (design Sect. 10)."""

    bands: dict[str, BandParams] = field(default_factory=lambda: dict(_DEFAULT_BANDS))
    depth_to_band: dict[str, str] = field(
        default_factory=lambda: dict(_DEFAULT_DEPTH_TO_BAND)
    )

    # T0
    gap_fail_hours: float = 6.0

    # T2a climatology
    k_clim: float = 4.0
    clim_doy_window_days: int = 15
    clim_min_years: int = 3

    # T2b ancillary consistency (ERA5 / LSM). depth-string -> ancillary column.
    ancillary_col: dict[str, str] = field(default_factory=dict)
    k_ancil_suspect: float = 4.0
    k_ancil_fail: float = 6.0
    ancil_resid_window_hours: float = 240.0  # 10-day rolling bias removal

    # T5 persistence
    repeat_suspect_hours: float = 6.0
    repeat_fail_hours: float = 12.0
    persist_win_hours: float = 24.0
    var_min_suspect: float = 0.05
    var_min_fail: float = 0.02
    t_freeze: float = 0.0
    delta_pc: float = 0.5  # -> 1.0 for saline / fine-textured soils

    # T6 profile consistency
    amp_gate: float = 4.0
    amp_tol: float = 0.15
    grad_suspect: dict[str, float] = field(
        default_factory=lambda: {"0-5": 1.5, "5-10": 1.0}
    )
    grad_fail: dict[str, float] = field(
        default_factory=lambda: {"0-5": 3.0, "5-10": 2.0}
    )


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def _dt_hours(index: pd.DatetimeIndex) -> float:
    """Native sampling interval (median spacing) in hours."""
    if len(index) < 2:
        return 1.0
    # total_seconds() is resolution-independent (pandas may store ns or us).
    deltas = index.to_series().diff().dt.total_seconds().to_numpy() / 3600.0
    deltas = deltas[np.isfinite(deltas) & (deltas > 0)]
    med = float(np.median(deltas)) if deltas.size else 1.0
    return med if med > 0 else 1.0


def _hours_to_samples(
    hours: float, dt_hours: float, *, odd: bool = False, floor: int = 1
) -> int:
    """Convert a physical window length to a sample count. ``floor`` guarantees a
    usable window even on coarse (e.g. 3-hourly) sites where the physical length
    would otherwise collapse to a single sample."""
    n = max(floor, int(round(hours / dt_hours)))
    if odd and n % 2 == 0:
        n += 1
    return n


def _min_periods(n: int, want: int) -> int:
    """Rolling ``min_periods`` clamped to never exceed the window length."""
    return min(max(1, want), n)


def _segment_ids(index: pd.DatetimeIndex, gap_fail_hours: float) -> np.ndarray:
    """Contiguous-segment id per sample. A time gap >= ``gap_fail_hours`` starts a
    new segment so shape/temporal tests (T3-T6) never span it."""
    gap_h = index.to_series().diff().dt.total_seconds().to_numpy() / 3600.0
    new_seg = np.isnan(gap_h) | (gap_h >= gap_fail_hours)
    return np.cumsum(new_seg) - 1


def _band_of(depth: str, config: QCConfig) -> str:
    band = config.depth_to_band.get(depth)
    if band is None:
        raise KeyError(
            f"depth {depth!r} not in depth_to_band {sorted(config.depth_to_band)}"
        )
    return band


def _mad(x: np.ndarray) -> float:
    """Median absolute deviation, NaN-aware."""
    x = x[~np.isnan(x)]
    if x.size == 0:
        return np.nan
    return float(np.median(np.abs(x - np.median(x))))


# --------------------------------------------------------------------------- #
# T0 - Structural
# --------------------------------------------------------------------------- #
def check_structural(values: pd.Series) -> np.ndarray:
    """MISSING for null, FAIL for non-finite (inf), PASS otherwise."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, PASS, dtype=np.int8)
    flag[~np.isfinite(v)] = FAIL
    flag[np.isnan(v)] = MISSING
    return flag


# --------------------------------------------------------------------------- #
# T1 - Gross / sensor range (static)
# --------------------------------------------------------------------------- #
def check_gross_range(values: pd.Series, band: BandParams) -> np.ndarray:
    """FAIL outside static physical bounds; MISSING carried through."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, PASS, dtype=np.int8)
    out = (v < band.gross_min) | (v > band.gross_max)
    flag[out] = FAIL
    flag[np.isnan(v)] = MISSING
    return flag


# --------------------------------------------------------------------------- #
# T2a - Climatology range (per-depth, per-DOY, multi-year)
# --------------------------------------------------------------------------- #
def check_climatology(values: pd.Series, config: QCConfig) -> np.ndarray:
    """Per-DOY envelope from this depth's own multi-year record. SUSPECT outside
    the empirical 0.5/99.5 percentile band or beyond ``K_CLIM`` sigma. Whole test
    is NOT_EVALUATED when the record is too short (< ``clim_min_years``)."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, PASS, dtype=np.int8)
    index = values.index

    n_years = index.year.nunique()
    if n_years < config.clim_min_years:
        flag[:] = NOT_EVALUATED
        flag[np.isnan(v)] = MISSING
        return flag

    doy = index.dayofyear.to_numpy()
    valid = ~np.isnan(v)

    # Per-DOY statistics, then smoothed over a circular +/- window.
    stats = _doy_envelope(doy[valid], v[valid], config)
    lo = stats["p0_5"][doy - 1]
    hi = stats["p99_5"][doy - 1]
    mean = stats["mean"][doy - 1]
    sd = stats["sd"][doy - 1]

    suspect = valid & (
        (v < lo) | (v > hi) | (np.abs(v - mean) > config.k_clim * sd)
    )
    flag[suspect] = SUSPECT
    # DOYs with no climatology support cannot be judged.
    flag[valid & np.isnan(mean)] = NOT_EVALUATED
    flag[np.isnan(v)] = MISSING
    return flag


def _doy_envelope(doy: np.ndarray, vals: np.ndarray, config: QCConfig) -> dict:
    """Circular +/- ``clim_doy_window_days`` percentile/mean/sd envelope, length 366."""
    w = config.clim_doy_window_days
    p0_5 = np.full(366, np.nan)
    p99_5 = np.full(366, np.nan)
    mean = np.full(366, np.nan)
    sd = np.full(366, np.nan)
    # Pre-bin values by DOY for cheap window aggregation.
    by_doy: dict[int, list] = {}
    for d, val in zip(doy, vals):
        by_doy.setdefault(int(d), []).append(val)

    for center in range(1, 367):
        offsets = [((center - 1 + k) % 366) + 1 for k in range(-w, w + 1)]
        pool = np.concatenate(
            [np.asarray(by_doy[o]) for o in offsets if o in by_doy]
        ) if any(o in by_doy for o in offsets) else np.empty(0)
        if pool.size >= 2:
            p0_5[center - 1] = np.percentile(pool, 0.5)
            p99_5[center - 1] = np.percentile(pool, 99.5)
            mean[center - 1] = pool.mean()
            sd[center - 1] = pool.std(ddof=1)
    return {"p0_5": p0_5, "p99_5": p99_5, "mean": mean, "sd": sd}


# --------------------------------------------------------------------------- #
# T2b - Ancillary consistency (ERA5 / land-surface model)
# --------------------------------------------------------------------------- #
def check_ancillary(
    values: pd.Series, ancillary: pd.Series, dt_hours: float, config: QCConfig
) -> np.ndarray:
    """Flag divergence of the sensor from an independent model soil temperature.

    The ancillary series is a *reference*, not truth (coarse grid cell, depth
    mismatch, systematic offset), so a slowly-varying bias is removed with a
    rolling-median of the residual and only a sudden robust departure is flagged.
    Capped at SUSPECT -- a model disagreement never rejects data on its own."""
    v = values.to_numpy(dtype=float)
    a = ancillary.to_numpy(dtype=float)
    flag = np.full(v.shape, PASS, dtype=np.int8)

    resid = v - a
    valid = np.isfinite(resid)
    if valid.sum() < 2:
        flag[:] = NOT_EVALUATED
        flag[np.isnan(v)] = MISSING
        return flag

    n = _hours_to_samples(config.ancil_resid_window_hours, dt_hours, odd=True, floor=3)
    mp = _min_periods(n, max(2, n // 4))
    r = pd.Series(resid)
    med = r.rolling(n, center=True, min_periods=mp).median()
    # Robust rolling scale of the de-biased residual.
    dev = (r - med).abs()
    scale = 1.4826 * dev.rolling(n, center=True, min_periods=mp).median()
    scale = scale.to_numpy()
    z = np.divide(
        dev.to_numpy(), scale, out=np.zeros_like(scale), where=scale > 0
    )

    ok = valid & np.isfinite(z)
    flag[ok & (z > config.k_ancil_suspect)] = SUSPECT  # SUSPECT ceiling
    flag[~valid] = NOT_EVALUATED  # ancillary missing where sensor present
    flag[np.isnan(v)] = MISSING
    return flag


# --------------------------------------------------------------------------- #
# T3 - Step / rate-of-change (segment-aware)
# --------------------------------------------------------------------------- #
def check_step(
    values: pd.Series, segment: np.ndarray, band: BandParams
) -> np.ndarray:
    """Per-hour first difference; segment boundaries emit NOT_EVALUATED."""
    v = values.to_numpy(dtype=float)
    index = values.index
    flag = np.full(v.shape, PASS, dtype=np.int8)

    dt_h = index.to_series().diff().dt.total_seconds().to_numpy() / 3600.0
    dv = np.abs(np.diff(v, prepend=np.nan))
    roc = np.where(dt_h > 0, dv / dt_h, np.nan)

    # First sample of each segment has no valid predecessor within the segment.
    seg_start = np.concatenate(([True], segment[1:] != segment[:-1]))
    roc[seg_start] = np.nan

    flag[roc > band.roc_suspect] = SUSPECT
    flag[roc > band.roc_fail] = FAIL
    flag[np.isnan(roc)] = NOT_EVALUATED
    flag[np.isnan(v)] = MISSING
    return flag


# --------------------------------------------------------------------------- #
# T4 - Spike (Hampel / MAD, noise-floored)
# --------------------------------------------------------------------------- #
def check_spike(
    values: pd.Series, segment: np.ndarray, dt_hours: float, band: BandParams
) -> np.ndarray:
    """Robust local-outlier test. The scale estimate is floored at
    ``SIGMA_NOISE`` so an isothermal plateau does not manufacture spikes and
    shred the zero-curtain (design Sect. 6 -- mandatory)."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, NOT_EVALUATED, dtype=np.int8)
    n = _hours_to_samples(band.spike_win_hours, dt_hours, odd=True, floor=3)
    mp = _min_periods(n, 3)

    for seg in np.unique(segment):
        idx = np.where(segment == seg)[0]
        if idx.size == 0:
            continue
        s = pd.Series(v[idx])
        med = s.rolling(n, center=True, min_periods=mp).median()
        mad = s.rolling(n, center=True, min_periods=mp).apply(_mad, raw=True)
        scale = np.maximum(1.4826 * mad.to_numpy(), band.sigma_noise)
        dev = np.abs(v[idx] - med.to_numpy())
        sub = flag[idx]
        evaluable = np.isfinite(dev) & np.isfinite(scale)
        sub[evaluable] = PASS
        sub[evaluable & (dev > band.n_spike_suspect * scale)] = SUSPECT
        sub[evaluable & (dev > band.n_spike_fail * scale)] = FAIL
        flag[idx] = sub

    flag[np.isnan(v)] = MISSING
    return flag


# --------------------------------------------------------------------------- #
# T5 - Persistence (split: T5a exact-repeat, T5b low-variance FT-gated)
# --------------------------------------------------------------------------- #
def check_persist_repeat(
    values: pd.Series, segment: np.ndarray, dt_hours: float, config: QCConfig
) -> np.ndarray:
    """T5a: consecutive bit-identical values = stuck logger. Always on,
    value-independent -- catches true faults even inside a 0 degC plateau, where a
    genuine curtain still carries small real fluctuation + digitization noise."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, PASS, dtype=np.int8)
    n_suspect = _hours_to_samples(config.repeat_suspect_hours, dt_hours)
    n_fail = _hours_to_samples(config.repeat_fail_hours, dt_hours)

    for seg in np.unique(segment):
        idx = np.where(segment == seg)[0]
        vv = v[idx]
        # Run id increments whenever the value changes or a NaN appears.
        change = np.concatenate(([True], (vv[1:] != vv[:-1]) | np.isnan(vv[1:])))
        run_id = np.cumsum(change)
        run_len = pd.Series(np.ones(vv.size)).groupby(run_id).cumcount() + 1
        # Total length of each run, broadcast back to every member.
        totals = pd.Series(run_len).groupby(run_id).transform("max").to_numpy()
        sub = flag[idx]
        finite = np.isfinite(vv)
        sub[finite & (totals >= n_suspect)] = SUSPECT
        sub[finite & (totals >= n_fail)] = FAIL
        flag[idx] = sub

    flag[np.isnan(v)] = MISSING
    return flag


def check_persist_var(
    values: pd.Series, segment: np.ndarray, dt_hours: float, config: QCConfig
) -> np.ndarray:
    """T5b: low rolling variance = stuck sensor -- but suppressed (NOT_EVALUATED)
    inside the phase-change window ``|T - T_FREEZE| < DELTA_PC``, where a flat run
    is the physical zero-curtain, not a fault."""
    v = values.to_numpy(dtype=float)
    flag = np.full(v.shape, NOT_EVALUATED, dtype=np.int8)
    n = _hours_to_samples(config.persist_win_hours, dt_hours, odd=True, floor=3)
    mp = _min_periods(n, 3)

    in_curtain = np.abs(v - config.t_freeze) < config.delta_pc

    for seg in np.unique(segment):
        idx = np.where(segment == seg)[0]
        s = pd.Series(v[idx])
        sd = s.rolling(n, center=True, min_periods=mp).std().to_numpy()
        sub = flag[idx]
        curtain = in_curtain[idx]
        evaluable = np.isfinite(sd) & ~curtain  # gated off inside the curtain
        sub[evaluable] = PASS
        sub[evaluable & (sd < config.var_min_suspect)] = SUSPECT
        sub[evaluable & (sd < config.var_min_fail)] = FAIL
        flag[idx] = sub

    flag[np.isnan(v)] = MISSING
    return flag


def combine_persistence(repeat: np.ndarray, var: np.ndarray) -> np.ndarray:
    """T5 combined = max(T5a, T5b). MISSING dominates."""
    out = np.maximum(repeat, var).astype(np.int8)
    out[(repeat == MISSING) | (var == MISSING)] = MISSING
    return out


# --------------------------------------------------------------------------- #
# T6 - Profile consistency (>= 2 depths; forcing-gated). Opportunistic.
# --------------------------------------------------------------------------- #
def check_profile_gradient(
    depth_values: dict[str, pd.Series], config: QCConfig
) -> dict[str, np.ndarray]:
    """T6c only (resolution-independent): vertical gradient bounds per adjacent
    pair. T6a (amplitude damping) / T6b (phase lag) need forcing/resolution gating
    and are left for a follow-up; the gradient check is the always-safe subset.

    Returns a per-depth flag aligned to each depth's own index. A pair violation
    flags *both* members (which sensor is wrong is decided by other tests)."""
    present = list(depth_values)
    flags = {
        d: np.full(len(depth_values[d]), NOT_EVALUATED, dtype=np.int8)
        for d in present
    }
    if len(present) < 2:
        return flags

    # Map present depths to nominal cm for pairing.
    order = sorted(present, key=_depth_cm)
    for upper, lower in zip(order, order[1:]):
        cu, cl = _depth_cm(upper), _depth_cm(lower)
        dz = abs(cl - cu)
        pair = f"{int(min(cu, cl))}-{int(max(cu, cl))}"
        g_susp = config.grad_suspect.get(pair)
        g_fail = config.grad_fail.get(pair)
        if g_susp is None or g_fail is None or dz == 0:
            continue

        joined = pd.concat(
            [depth_values[upper].rename("u"), depth_values[lower].rename("l")],
            axis=1, join="inner",
        )
        grad = (joined["u"] - joined["l"]).abs() / dz
        f = np.full(len(joined), PASS, dtype=np.int8)
        f[grad.to_numpy() > g_susp] = SUSPECT
        f[grad.to_numpy() > g_fail] = FAIL
        f[~np.isfinite(grad.to_numpy())] = NOT_EVALUATED

        for d in (upper, lower):
            pos = depth_values[d].index.get_indexer(joined.index)
            keep = pos >= 0
            worse = np.maximum(flags[d][pos[keep]], f[keep])
            # NOT_EVALUATED default should yield to a real PASS/higher result.
            base = flags[d][pos[keep]]
            base_eval = np.where(base == NOT_EVALUATED, f[keep], worse)
            flags[d][pos[keep]] = base_eval
    return flags


def _depth_cm(depth: str) -> float:
    """Parse ``"5cm"`` / ``"10cm"`` -> 5.0 / 10.0."""
    return float("".join(c for c in depth if (c.isdigit() or c == ".")) or "nan")


# --------------------------------------------------------------------------- #
# T7 - Aggregation
# --------------------------------------------------------------------------- #
def aggregate(flag_columns: dict[str, np.ndarray]) -> np.ndarray:
    """Summary flag = worst-of (max) across per-test flags. 9 dominates when the
    value is missing; the ordinal ordering makes plain ``max`` correct."""
    stacked = np.vstack(list(flag_columns.values()))
    summary = stacked.max(axis=0).astype(np.int8)
    if np.any(stacked == MISSING):
        summary[np.any(stacked == MISSING, axis=0)] = MISSING
    return summary


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def _discover_depths(df: pd.DataFrame, config: QCConfig) -> dict[str, str]:
    """Map depth-string -> column name for every ``soil_temp_<depth>`` present
    that also has a band mapping. E.g. ``soil_temp_5cm`` -> ``{"5cm": ...}``."""
    found: dict[str, str] = {}
    for col in df.columns:
        if col.startswith("soil_temp_"):
            depth = col[len("soil_temp_"):]
            if depth in config.depth_to_band:
                found[depth] = col
    return found


def apply_qc(
    df: pd.DataFrame, config: QCConfig | None = None
) -> tuple[pd.DataFrame, dict]:
    """Run the full T0-T7 QC pipeline over every discovered depth.

    Input: a DataFrame with a ``DatetimeIndex`` and one or more
    ``soil_temp_<depth>`` columns (e.g. ``soil_temp_5cm``). Optional ancillary
    model-temperature columns are wired per depth via ``config.ancillary_col``.

    Returns ``(out, report)`` where ``out`` is the input frame with per-test flag
    columns ``qf_<test>_<depth>`` and a roll-up ``qf_summary_<depth>`` added
    (raw values untouched), and ``report`` holds run diagnostics in ``.attrs``.
    """
    if config is None:
        config = QCConfig()
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError("apply_qc requires a DatetimeIndex.")

    out = df.copy()
    index = out.index
    dt_hours = _dt_hours(index)
    segment = _segment_ids(index, config.gap_fail_hours)

    depths = _discover_depths(out, config)
    if not depths:
        raise ValueError(
            "No soil_temp_<depth> columns found matching config.depth_to_band."
        )

    depth_values = {d: out[col].astype(float) for d, col in depths.items()}
    profile_grad = check_profile_gradient(depth_values, config)

    report: dict = {"dt_hours": dt_hours, "n_segments": int(segment.max()) + 1,
                    "depths": {}}

    for depth, col in depths.items():
        band = config.bands[_band_of(depth, config)]
        vals = depth_values[depth]

        per_test: dict[str, np.ndarray] = {}
        per_test["qf_struct"] = check_structural(vals)
        per_test["qf_range"] = check_gross_range(vals, band)
        per_test["qf_clim"] = check_climatology(vals, config)
        per_test["qf_step"] = check_step(vals, segment, band)
        per_test["qf_spike"] = check_spike(vals, segment, dt_hours, band)

        repeat = check_persist_repeat(vals, segment, dt_hours, config)
        var = check_persist_var(vals, segment, dt_hours, config)
        per_test["qf_persist_repeat"] = repeat
        per_test["qf_persist_var"] = var
        per_test["qf_persist"] = combine_persistence(repeat, var)

        # T2b ancillary consistency (optional, only if a column is configured).
        anc_col = config.ancillary_col.get(depth)
        if anc_col and anc_col in out.columns:
            per_test["qf_ancil"] = check_ancillary(
                vals, out[anc_col].astype(float), dt_hours, config
            )

        # T6 profile gradient (only meaningful with >=2 depths).
        if len(depths) >= 2:
            per_test["qf_profile"] = profile_grad[depth]

        # Write per-test columns.
        for name, arr in per_test.items():
            out[f"{name}_{depth}"] = arr

        # Summary excludes the split sub-flags (kept for provenance only).
        summary_inputs = {
            k: v for k, v in per_test.items()
            if k not in ("qf_persist_repeat", "qf_persist_var")
        }
        out[f"qf_summary_{depth}"] = aggregate(summary_inputs)
        report["depths"][depth] = {
            "band": _band_of(depth, config),
            "tests": sorted(summary_inputs),
        }

    out.attrs["qc"] = report
    return out, report
