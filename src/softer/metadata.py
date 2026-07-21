"""Per-site metadata loader.

Reads the site-metadata CSV that holds per-site facts used by the pipeline —
most importantly ``sigma_t`` (the instrument-specific soil-temperature
uncertainty that sets the freeze/thaw deadband). Anything constant across sites
lives in the YAML config instead (:mod:`softer.config`).

Expected CSV columns (extra columns are preserved and ignored here):
    site_id, sensor, sigma_t

If a site row omits ``sigma_t``, it falls back to the sensor's default from the
:mod:`softer.sensors` registry.
"""

from __future__ import annotations

import math

import pandas as pd

from .sensors import get_sigma_t as _sensor_sigma_t


def load_site_metadata(path: str) -> pd.DataFrame:
    """Load the site-metadata CSV, indexed by ``site_id``."""
    meta = pd.read_csv(path)
    if "site_id" not in meta.columns:
        raise KeyError("Site-metadata CSV must have a 'site_id' column.")
    meta = meta.set_index("site_id")
    if meta.index.duplicated().any():
        dupes = meta.index[meta.index.duplicated()].unique().tolist()
        raise ValueError(f"Duplicate site_id(s) in metadata: {dupes}.")
    return meta


def get_sigma_t(meta: pd.DataFrame, site_id: str) -> float:
    """Return sigma_t for a site: explicit value if present, else the sensor default."""
    if site_id not in meta.index:
        raise KeyError(f"Site {site_id!r} not found in metadata.")
    row = meta.loc[site_id]

    sigma = row.get("sigma_t") if hasattr(row, "get") else None
    if sigma is not None and not (isinstance(sigma, float) and math.isnan(sigma)):
        return float(sigma)

    sensor = row.get("sensor") if hasattr(row, "get") else None
    if sensor is not None and not (isinstance(sensor, float) and math.isnan(sensor)):
        return _sensor_sigma_t(str(sensor))

    raise ValueError(
        f"Site {site_id!r} has neither a sigma_t value nor a known sensor to fall back on."
    )
