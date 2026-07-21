"""Tests for config and site-metadata loaders."""

import math

import pandas as pd
import pytest

from softer.config import load_config
from softer.metadata import get_sigma_t, load_site_metadata


def test_load_example_config():
    cfg = load_config("configs/config.example.yaml")
    assert cfg.cycle_detection.year_start == (8, 1)
    assert cfg.cycle_detection.min_duration == "2D"
    assert cfg.usability.warm_edge == 2.0
    assert cfg.usability.cold_edge == -2.0
    assert cfg.usability.min_reach == -1.0
    assert cfg.usability.max_temp_gap == 0.5
    assert cfg.coverage.critical_window.start == (10, 1)
    assert cfg.coverage.critical_window.end == (3, 1)
    assert cfg.coverage.min_coverage == 0.80


def test_config_defaults_and_unknown_key(tmp_path):
    # Empty config -> all defaults.
    p = tmp_path / "empty.yaml"
    p.write_text("")
    cfg = load_config(str(p))
    assert cfg.usability.max_temp_gap == 0.5

    # Unknown key -> error.
    bad = tmp_path / "bad.yaml"
    bad.write_text("usability:\n  bogus: 1\n")
    with pytest.raises(KeyError):
        load_config(str(bad))


def test_metadata_sigma_and_sensor_fallback():
    meta = load_site_metadata("configs/site_metadata.example.csv")
    # Explicit sigma_t.
    assert get_sigma_t(meta, "FM403") == 0.375
    # Missing sigma_t falls back to the sensor default (TEROS12).
    gr = get_sigma_t(meta, "GR01")
    assert isinstance(gr, float) and not math.isnan(gr)


def test_metadata_missing_site():
    meta = load_site_metadata("configs/site_metadata.example.csv")
    with pytest.raises(KeyError):
        get_sigma_t(meta, "NOPE")
