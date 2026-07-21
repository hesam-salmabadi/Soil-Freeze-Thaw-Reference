"""Tests for freeze/thaw cycle labeling (softer.preprocess.cycles)."""

import numpy as np
import pandas as pd

from softer.preprocess.cycles import (
    FREEZING,
    NO_CYCLE,
    THAWING,
    label_freeze_thaw_cycles,
)


def _series_to_frame(temps):
    idx = pd.date_range("2021-09-01", periods=len(temps), freq="1h")
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

    # Both phases present, exactly two cycles (one freeze, one thaw).
    assert set(out["cycle_phase"].unique()) >= {FREEZING, THAWING}
    assert out.loc[out["cycle_id"] >= 0, "cycle_id"].nunique() == 2

    # Pre-peak warm rise is outside any cycle.
    assert out["cycle_phase"].iloc[0] == NO_CYCLE

    # Descent into the global minimum is freezing; ascent out of it is thawing.
    tmin_pos = int(np.argmin(out["Temp"].values))
    assert out["cycle_phase"].iloc[tmin_pos - 1] == FREEZING
    assert out["cycle_phase"].iloc[tmin_pos + 1] == THAWING

    assert out.attrs["never_frozen"] is False


def test_near_zero_wiggles_do_not_spawn_cycles():
    # A single legitimate crossing, but jittery around zero within the deadband.
    temps = (
        [1, 2, 3, 4, 5]                       # warm peak
        + [2, 1.5, 0.3, -0.3, 0.3, -0.3, -1.5]  # jitter within +/-0.5 around the crossing
        + [-3, -5, -8]                        # clear trough
        + [-4, -1, 2, 4]                      # thaw back up
    )
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    # Despite the near-zero jitter, only two cycles should be detected.
    assert out.loc[out["cycle_id"] >= 0, "cycle_id"].nunique() == 2


def test_never_frozen_series():
    temps = [2, 3, 4, 5, 4, 3, 2, 3, 4]  # always well above zero
    df = _series_to_frame(temps)

    out = label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5)

    assert (out["cycle_phase"] == NO_CYCLE).all()
    assert (out["cycle_id"] == NO_CYCLE).all()
    assert out.attrs["never_frozen"] is True


def test_passthrough_columns_and_time_col():
    temps = [5, 2, -2, -6, -2, 2, 5]
    idx = pd.date_range("2021-10-01", periods=len(temps), freq="1h")
    df = pd.DataFrame(
        {
            "timestamp": idx,
            "Temp": temps,
            "EDC": np.linspace(20, 5, len(temps)),  # moisture rides along untouched
        }
    )

    out = label_freeze_thaw_cycles(df, temp_col="Temp", time_col="timestamp", deadband_c=0.5)

    # Moisture column preserved, new label columns added.
    assert "EDC" in out.columns
    assert {"cycle_id", "cycle_phase"}.issubset(out.columns)
    assert len(out) == len(df)


def test_deadband_from_sensor():
    temps = [5, 2, -2, -6, -2, 2, 5]
    df = _series_to_frame(temps)

    # Should resolve the deadband from the sensor registry without error.
    out = label_freeze_thaw_cycles(df, temp_col="Temp", sensor="TEROS12")
    assert {"cycle_id", "cycle_phase"}.issubset(out.columns)
