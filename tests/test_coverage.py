"""Tests for the usability / temperature-coverage gate (softer.preprocess.coverage)."""

import numpy as np
import pandas as pd

from softer.preprocess.cycles import label_freeze_thaw_cycles
from softer.preprocess.coverage import assess_usability, write_report


def _label(temps, start="2021-10-01", freq="1h", **kw):
    idx = pd.date_range(start, periods=len(temps), freq=freq)
    df = pd.DataFrame({"Temp": temps}, index=idx)
    return label_freeze_thaw_cycles(df, temp_col="Temp", deadband_c=0.5, max_gap=None,
                                    min_duration=None, **kw)


def test_well_sampled_freezing_cycle_is_usable():
    # Dense descent from +2 through 0 to -3, then back up.
    down = list(np.arange(2.0, -3.01, -0.2))
    up = list(np.arange(-3.0, 2.01, 0.2))
    labeled = _label(down + up)
    out, report = assess_usability(labeled, temp_col="Temp", warm_edge=2.0,
                                   cold_edge=-2.0, min_reach=-1.0, max_temp_gap=0.5)
    assert (report["usable"]).any()
    # The freezing rows of a usable cycle are marked True.
    assert out.loc[out["cycle_phase"] == 1, "usable"].any()


def test_temperature_hole_makes_cycle_unusable():
    # Jump straight from +2 to -3 (a >0.5C hole through the transition), then up.
    down = [2.0, -3.0]
    up = list(np.arange(-3.0, 2.01, 0.2))
    labeled = _label(down + up)
    out, report = assess_usability(labeled, temp_col="Temp", warm_edge=2.0,
                                   cold_edge=-2.0, min_reach=-1.0, max_temp_gap=0.5)
    freezing = report  # freezing cycles only
    # The freezing cycle with the big transition hole is not usable.
    assert not bool(freezing.loc[freezing["cycle_id"] == freezing["cycle_id"].min(), "usable"].iloc[0])


def test_shallow_cycle_flagged_and_not_usable():
    # Only reaches -0.6 (warmer than min_reach -1) -> shallow -> not usable.
    down = list(np.arange(2.0, -0.61, -0.2))
    up = list(np.arange(-0.6, 2.01, 0.2))
    labeled = _label(down + up)
    out, report = assess_usability(labeled, temp_col="Temp", warm_edge=2.0,
                                   cold_edge=-2.0, min_reach=-1.0, max_temp_gap=0.5)
    assert bool(report["shallow"].any())
    assert not bool(report["usable"].any())


def test_report_has_expected_columns_and_writes(tmp_path):
    down = list(np.arange(2.0, -3.01, -0.2))
    up = list(np.arange(-3.0, 2.01, 0.2))
    labeled = _label(down + up)
    out, report = assess_usability(labeled, temp_col="Temp", site_id="FM403")
    for col in ("site_id", "freezing_year", "cycle_id", "usable", "shallow",
                "max_temp_hole", "time_coverage", "longest_gap_days", "reason"):
        assert col in report.columns

    path = tmp_path / "usability.csv"
    write_report(report, str(path))
    back = pd.read_csv(path)
    assert len(back) == len(report)
    assert back["site_id"].iloc[0] == "FM403"
