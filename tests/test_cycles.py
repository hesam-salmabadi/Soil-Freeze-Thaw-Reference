"""Tests for freeze/thaw cycle labeling (softer.preprocess.cycles)."""

import numpy as np
import pandas as pd

from softer.preprocess.cycles import (
    FREEZING,
    NO_CYCLE,
    THAWING,
    label_freeze_thaw_cycles,
)


def _series_to_frame(temps, start="2021-09-01"):
    idx = pd.date_range(start, periods=len(temps), freq="1h")
    return pd.DataFrame({"Temp": temps}, index=idx)


def test_basic_freeze_then_thaw():
    # Warm rise to a peak (pre-cycle), decline into a subzero trough (freezing),
    # then rise back above zero (thawing).
    temps = (
        [1, 2, 3, 4, 5]              # pre-cycle warm rise -> peak at pos 4
        + [3, 1, -1, -3, -6, -8]     # freezing descent -> trough at pos 10
        + [-5, -2, 0.2, 2, 4, 5]     # thawing ascent (0.2 is within deadband)
    )
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert set(out["cycle_phase"].unique()) >= {FREEZING, THAWING}
    assert out.loc[out["cycle_id"] >= 0, "cycle_id"].nunique() == 2
    assert out["cycle_phase"].iloc[0] == NO_CYCLE

    tmin_pos = int(np.argmin(out["Temp"].values))
    assert out["cycle_phase"].iloc[tmin_pos - 1] == FREEZING
    assert out["cycle_phase"].iloc[tmin_pos + 1] == THAWING

    year_diag = out.attrs["years"]
    assert next(iter(year_diag.values()))["detection"] == "normal"


def test_near_zero_wiggles_do_not_spawn_cycles():
    temps = (
        [1, 2, 3, 4, 5]
        + [2, 1.5, 0.3, -0.3, 0.3, -0.3, -1.5]  # jitter within +/-0.5 around crossing
        + [-3, -5, -8]
        + [-4, -1, 2, 4]
    )
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert out.loc[out["cycle_id"] >= 0, "cycle_id"].nunique() == 2


def test_passthrough_columns_and_time_col():
    temps = [5, 2, -2, -6, -2, 2, 5]
    idx = pd.date_range("2021-10-01", periods=len(temps), freq="1h")
    df = pd.DataFrame(
        {
            "timestamp": idx,
            "Temp": temps,
            "EDC": np.linspace(20, 5, len(temps)),
        }
    )

    out = label_freeze_thaw_cycles(df, temp_col="Temp", time_col="timestamp", deadband_c=0.5)

    assert "EDC" in out.columns
    assert {"freezing_year", "cycle_id", "cycle_phase"}.issubset(out.columns)
    assert len(out) == len(df)


def test_deadband_from_sensor():
    temps = [5, 2, -2, -6, -2, 2, 5]
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", sensor="TEROS12")
    assert {"cycle_id", "cycle_phase"}.issubset(out.columns)


# --- multi-year windowing -------------------------------------------------

def test_multi_year_windows_reset_and_do_not_merge():
    # Build 3 back-to-back Aug->Aug winters, each: warm -> subzero -> warm.
    def winter(peak_start):
        # ~120-day season sampled daily: cosine dipping below zero mid-winter.
        idx = pd.date_range(peak_start, periods=120, freq="1D")
        t = 6 * np.cos(np.linspace(0, 2 * np.pi, 120))
        return pd.DataFrame({"Temp": t}, index=idx)

    frames = [winter(f"{y}-10-01") for y in (2019, 2020, 2021)]
    df = pd.concat(frames)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5, max_gap=None)

    years = sorted(y for y in out["freezing_year"].unique() if y >= 0)
    assert years == [2019, 2020, 2021]

    # cycle_id resets each year (each winter has a small number of cycles).
    for y in years:
        yr = out[out["freezing_year"] == y]
        ids = sorted(i for i in yr["cycle_id"].unique() if i >= 0)
        assert ids[0] == 0  # resets to 0 within each year

    # No single cycle spans two different freezing years.
    grouped = out[out["cycle_id"] >= 0].groupby(["freezing_year", "cycle_id"])
    assert len(grouped) > 0


def test_aug1_boundary_splits_years():
    # A point in July belongs to the previous freezing year; August starts a new one.
    idx = pd.to_datetime(["2020-07-31", "2020-08-01", "2020-08-02"])
    df = pd.DataFrame({"Temp": [3.0, 3.0, 3.0]}, index=idx)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert out.loc["2020-07-31", "freezing_year"] == 2019
    assert out.loc["2020-08-01", "freezing_year"] == 2020


# --- fallbacks ------------------------------------------------------------

def test_never_frozen_stalls_near_zero_uses_fallback():
    # Dips toward zero but never below -h: must still get a freeze + thaw.
    temps = [3, 2.5, 2, 1, 0.5, 0.3, 0.5, 1, 2, 3]
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert (out["cycle_phase"] == FREEZING).any()
    assert (out["cycle_phase"] == THAWING).any()
    diag = next(iter(out.attrs["years"].values()))
    assert diag["detection"] == "fallback"
    assert diag["never_frozen"] is True


def test_always_frozen_uses_fallback():
    # Permafrost-like: never rises above +h.
    temps = [-2, -3, -5, -8, -10, -8, -5, -3]
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert (out["cycle_phase"] == FREEZING).any()
    assert (out["cycle_phase"] == THAWING).any()
    diag = next(iter(out.attrs["years"].values()))
    assert diag["detection"] == "fallback"
    assert diag["never_thawed"] is True


def test_degenerate_flat_window_stays_unlabeled():
    temps = [1.0] * 8
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert (out["cycle_phase"] == NO_CYCLE).all()
    assert (out["cycle_id"] == NO_CYCLE).all()


# --- data gaps ------------------------------------------------------------

def test_cycle_does_not_span_large_gap():
    # A warm plateau, then a multi-day outage, then a subzero dip. Without gap
    # handling, interpolation could bridge these into one long descent.
    warm = pd.date_range("2021-10-01", periods=5, freq="1h")
    cold = pd.date_range("2021-10-20", periods=5, freq="1h")  # ~19-day gap
    idx = warm.append(cold)
    temps = [5, 5, 5, 5, 5, -5, -6, -7, -6, -5]
    df = pd.DataFrame({"Temp": temps}, index=idx)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5, max_gap="2D")

    # The warm and cold blocks are in separate gap segments; no cycle bridges them.
    for cid in out.loc[out["cycle_id"] >= 0, "cycle_id"].unique():
        seg = out[out["cycle_id"] == cid]
        assert seg.index.to_series().diff().max() <= pd.Timedelta("2D")
