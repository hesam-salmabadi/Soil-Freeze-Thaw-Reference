"""Tests for soil-temperature QA/QC (softer.preprocess.quality_control)."""

import numpy as np
import pandas as pd

from softer.preprocess.quality_control import (
    FAIL,
    MISSING,
    NOT_EVALUATED,
    PASS,
    SUSPECT,
    QCConfig,
    aggregate,
    apply_qc,
    check_ancillary,
    check_climatology,
    check_gross_range,
    check_persist_repeat,
    check_persist_var,
    check_profile_gradient,
    check_spike,
    check_step,
    check_structural,
    combine_persistence,
)


def _series(values, start="2021-01-01", freq="1h"):
    idx = pd.date_range(start, periods=len(values), freq=freq)
    return pd.Series(np.asarray(values, dtype=float), index=idx)


def _seg(series):
    from softer.preprocess.quality_control import _segment_ids

    return _segment_ids(series.index, 6.0)


# --------------------------------------------------------------------------- #
# T0 - Structural
# --------------------------------------------------------------------------- #
def test_structural_missing_and_nonfinite():
    s = _series([1.0, np.nan, np.inf, -np.inf, 5.0])
    flag = check_structural(s)
    assert flag[0] == PASS
    assert flag[1] == MISSING
    assert flag[2] == FAIL   # +inf is non-finite but present
    assert flag[3] == FAIL
    assert flag[4] == PASS


# --------------------------------------------------------------------------- #
# T1 - Gross range
# --------------------------------------------------------------------------- #
def test_gross_range_bounds():
    band = QCConfig().bands["SHAL"]  # [-40, 55]
    s = _series([-41.0, 0.0, 55.1, 20.0])
    flag = check_gross_range(s, band)
    assert flag.tolist() == [FAIL, PASS, FAIL, PASS]


def test_surf_band_wider_than_shal():
    cfg = QCConfig()
    s = _series([65.0])  # inside SURF (+70) but above SHAL (+55)
    assert check_gross_range(s, cfg.bands["SURF"])[0] == PASS
    assert check_gross_range(s, cfg.bands["SHAL"])[0] == FAIL


# --------------------------------------------------------------------------- #
# T2a - Climatology
# --------------------------------------------------------------------------- #
def test_climatology_short_record_not_evaluated():
    # Single year -> below clim_min_years -> whole test NOT_EVALUATED.
    s = _series(np.zeros(24 * 30), freq="1h")
    flag = check_climatology(s, QCConfig())
    assert set(np.unique(flag)) == {NOT_EVALUATED}


def test_climatology_flags_seasonal_outlier():
    idx = pd.date_range("2018-01-01", periods=24 * 365 * 4, freq="1h")
    seas = 10 * np.cos((idx.dayofyear.to_numpy() - 30) / 365 * 2 * np.pi)
    rng = np.random.default_rng(0)
    v = seas + rng.normal(0, 0.5, len(idx))
    v[5000] = 40.0  # far outside the DOY envelope but inside gross range
    s = pd.Series(v, index=idx)
    flag = check_climatology(s, QCConfig())
    assert flag[5000] == SUSPECT
    assert PASS in np.unique(flag)


# --------------------------------------------------------------------------- #
# T2b - Ancillary consistency
# --------------------------------------------------------------------------- #
def test_ancillary_ignores_constant_bias():
    # Sensor tracks ERA5 with a fixed +2 degC offset -> no flags (bias removed).
    rng = np.random.default_rng(1)
    a = _series(5 * np.sin(np.arange(500) / 20) + rng.normal(0, 0.1, 500))
    v = a + 2.0
    flag = check_ancillary(v, a, 1.0, QCConfig())
    assert set(np.unique(flag)) <= {PASS}


def test_ancillary_flags_sudden_divergence_capped_at_suspect():
    rng = np.random.default_rng(2)
    a = _series(5 * np.sin(np.arange(500) / 20) + rng.normal(0, 0.1, 500))
    # Sensor tracks the model with its own realistic noise (nonzero residual
    # spread), then a single sample diverges hard.
    v = a + rng.normal(0, 0.15, 500)
    v.iloc[250] += 15.0
    flag = check_ancillary(v, a, 1.0, QCConfig())
    assert flag[250] == SUSPECT
    assert FAIL not in np.unique(flag)  # model disagreement never rejects


def test_ancillary_missing_where_model_absent():
    a = _series([1.0, np.nan, 3.0, 4.0, 5.0])
    v = _series([1.1, 2.0, 3.1, 4.1, 5.1])
    flag = check_ancillary(v, a, 1.0, QCConfig())
    assert flag[1] == NOT_EVALUATED


# --------------------------------------------------------------------------- #
# T3 - Step / rate-of-change
# --------------------------------------------------------------------------- #
def test_step_rate_of_change():
    band = QCConfig().bands["SHAL"]  # suspect 4, fail 6 degC/hr
    s = _series([0.0, 0.0, 4.5, 4.5, 11.0])  # jumps of 0, 4.5, 0, 6.5 per hour
    flag = check_step(s, _seg(s), band)
    assert flag[0] == NOT_EVALUATED       # first sample, no predecessor
    assert flag[2] == SUSPECT             # 4.5 > 4
    assert flag[4] == FAIL                # 6.5 > 6


def test_step_not_evaluated_across_gap():
    idx = pd.DatetimeIndex(
        ["2021-01-01 00:00", "2021-01-01 01:00", "2021-01-01 12:00"]
    )
    s = pd.Series([0.0, 0.0, 20.0], index=idx)  # 11h gap before the jump
    flag = check_step(s, _seg(s), QCConfig().bands["SHAL"])
    assert flag[2] == NOT_EVALUATED  # step must not span the gap


# --------------------------------------------------------------------------- #
# T4 - Spike
# --------------------------------------------------------------------------- #
def test_spike_detects_outlier():
    rng = np.random.default_rng(3)
    v = 1.0 + rng.normal(0, 0.05, 200)
    v[100] = 25.0
    s = _series(v)
    flag = check_spike(s, _seg(s), 1.0, QCConfig().bands["SHAL"])
    assert flag[100] == FAIL


def test_spike_noise_floor_protects_plateau():
    # A perfectly flat plateau with tiny noise below the floor must NOT spike.
    rng = np.random.default_rng(4)
    v = -0.3 + rng.normal(0, 0.02, 200)  # under SIGMA_NOISE(SHAL)=0.08
    s = _series(v)
    flag = check_spike(s, _seg(s), 1.0, QCConfig().bands["SHAL"])
    assert SUSPECT not in np.unique(flag)
    assert FAIL not in np.unique(flag)


# --------------------------------------------------------------------------- #
# T5 - Persistence (split)
# --------------------------------------------------------------------------- #
def test_persist_repeat_stuck_logger():
    v = np.concatenate([np.linspace(2, 1, 50), np.full(40, 0.5)])  # 40h identical
    s = _series(v)
    flag = check_persist_repeat(s, _seg(s), 1.0, QCConfig())
    assert flag[70] == FAIL  # deep inside the 40h stuck run (>= 12h)


def test_persist_var_gated_off_in_curtain():
    # Flat run pinned at the freezing point -> variance test suppressed.
    s = _series(np.full(60, -0.1))  # |T| < DELTA_PC(0.5)
    flag = check_persist_var(s, _seg(s), 1.0, QCConfig())
    assert set(np.unique(flag)) <= {NOT_EVALUATED, MISSING}


def test_persist_var_flags_stuck_away_from_freezing():
    # Flat run far from 0 degC -> low variance IS a stuck-sensor signal.
    s = _series(np.full(60, 10.0))
    flag = check_persist_var(s, _seg(s), 1.0, QCConfig())
    assert FAIL in np.unique(flag)


def test_combine_persistence_takes_worst():
    repeat = np.array([PASS, SUSPECT, PASS, MISSING], dtype=np.int8)
    var = np.array([PASS, PASS, FAIL, PASS], dtype=np.int8)
    out = combine_persistence(repeat, var)
    assert out.tolist() == [PASS, SUSPECT, FAIL, MISSING]


# --------------------------------------------------------------------------- #
# T6 - Profile gradient
# --------------------------------------------------------------------------- #
def test_profile_single_depth_not_evaluated():
    flags = check_profile_gradient({"5cm": _series([1.0, 2.0])}, QCConfig())
    assert set(np.unique(flags["5cm"])) == {NOT_EVALUATED}


def test_profile_gradient_flags_extreme_pair():
    # 0-5 cm pair: fail bound 3.0 degC/cm -> dz=5cm -> |dT|=20 -> 4 degC/cm.
    surf = _series([25.0, 25.0])
    shal = _series([5.0, 5.0])
    flags = check_profile_gradient({"0cm": surf, "5cm": shal}, QCConfig())
    assert flags["0cm"][0] == FAIL
    assert flags["5cm"][0] == FAIL


# --------------------------------------------------------------------------- #
# T7 - Aggregation
# --------------------------------------------------------------------------- #
def test_aggregate_worst_of_and_missing_dominates():
    cols = {
        "a": np.array([PASS, SUSPECT, PASS], dtype=np.int8),
        "b": np.array([PASS, PASS, MISSING], dtype=np.int8),
    }
    out = aggregate(cols)
    assert out.tolist() == [PASS, SUSPECT, MISSING]


# --------------------------------------------------------------------------- #
# Integration - apply_qc
# --------------------------------------------------------------------------- #
def _synthetic_frame():
    idx = pd.date_range("2021-01-01", periods=24 * 40, freq="1h")
    rng = np.random.default_rng(0)
    base = 2 * np.sin(np.arange(len(idx)) / 24 * 2 * np.pi * 0.1) + rng.normal(
        0, 0.1, len(idx)
    )
    t5 = base.copy()
    t10 = base * 0.6 + rng.normal(0, 0.05, len(idx))
    t5[100] = 25.0     # spike
    t5[200] = 999.0    # gross fail
    t5[300:340] = -0.30  # stuck run
    era5 = base + 0.5
    era5[100] = base[100]
    return pd.DataFrame(
        {"soil_temp_5cm": t5, "soil_temp_10cm": t10, "era5_stl1": era5}, index=idx
    )


def test_apply_qc_end_to_end():
    df = _synthetic_frame()
    cfg = QCConfig(ancillary_col={"5cm": "era5_stl1"})
    out, report = apply_qc(df, cfg)

    # Raw values untouched.
    assert (out["soil_temp_5cm"] == df["soil_temp_5cm"]).all()
    # Expected per-test columns exist for both depths.
    for depth in ("5cm", "10cm"):
        for test in ("struct", "range", "clim", "step", "spike", "persist",
                     "persist_repeat", "persist_var", "profile", "summary"):
            assert f"qf_{test}_{depth}" in out.columns
    # Ancillary column only where configured.
    assert "qf_ancil_5cm" in out.columns
    assert "qf_ancil_10cm" not in out.columns

    # Known injected faults are caught.
    assert out["qf_range_5cm"].iloc[200] == FAIL
    assert out["qf_spike_5cm"].iloc[100] == FAIL
    assert out["qf_persist_repeat_5cm"].iloc[300:340].max() == FAIL
    assert out["qf_summary_5cm"].iloc[200] == FAIL
    assert report["dt_hours"] == 1.0


def test_apply_qc_missing_values_flagged_missing():
    idx = pd.date_range("2021-01-01", periods=5, freq="1h")
    df = pd.DataFrame(
        {"soil_temp_5cm": [1.0, np.nan, 2.0, np.nan, 3.0]}, index=idx
    )
    out, _ = apply_qc(df, QCConfig())
    assert out["qf_summary_5cm"].tolist() == [
        NOT_EVALUATED, MISSING, NOT_EVALUATED, MISSING, NOT_EVALUATED
    ]


def test_apply_qc_three_hourly_with_gap():
    idx = pd.date_range("2021-01-01", periods=200, freq="3h")
    keep = np.ones(200, bool)
    keep[100:104] = False  # 12h hole -> new segment
    idx = idx[keep]
    rng = np.random.default_rng(1)
    v = 1.0 + rng.normal(0, 0.05, idx.size)
    out, report = apply_qc(pd.DataFrame({"soil_temp_5cm": v}, index=idx), QCConfig())
    assert report["dt_hours"] == 3.0
    assert report["n_segments"] == 2
    assert (out["qf_step_5cm"] == NOT_EVALUATED).sum() >= 2


def test_apply_qc_single_depth_has_no_profile():
    idx = pd.date_range("2021-01-01", periods=100, freq="1h")
    df = pd.DataFrame({"soil_temp_5cm": np.zeros(100)}, index=idx)
    out, _ = apply_qc(df, QCConfig())
    assert "qf_profile_5cm" not in out.columns


def test_apply_qc_requires_datetime_index():
    df = pd.DataFrame({"soil_temp_5cm": [1.0, 2.0]})
    try:
        apply_qc(df, QCConfig())
    except TypeError:
        return
    raise AssertionError("expected TypeError for non-DatetimeIndex input")


def test_apply_qc_no_depth_columns_raises():
    idx = pd.date_range("2021-01-01", periods=3, freq="1h")
    df = pd.DataFrame({"temperature": [1.0, 2.0, 3.0]}, index=idx)
    try:
        apply_qc(df, QCConfig())
    except ValueError:
        return
    raise AssertionError("expected ValueError when no soil_temp_<depth> columns")
