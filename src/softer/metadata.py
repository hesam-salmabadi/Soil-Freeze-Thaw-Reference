"""Per-site metadata loader.

Reads the site-metadata CSV that holds per-site facts used by the pipeline —
most importantly ``sigma_t`` (the instrument-specific soil-temperature
uncertainty that sets the freeze/thaw deadband). Anything constant across sites
lives in the YAML config instead (:mod:`softer.config`).

Expected CSV columns (extra columns are preserved and ignored here):
    site_id, sensor, sigma_t

``sigma_t`` is the single source of truth for the deadband. If a site row leaves
it blank, :func:`get_sigma_t` returns ``None`` and the caller falls back to the
configured ``default_deadband``.
"""

from __future__ import annotations

import math

import pandas as pd


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


def get_sigma_t(meta: pd.DataFrame, site_id: str) -> float | None:
    """Return the site's ``sigma_t`` from metadata, or ``None`` if it is blank.

    A ``None`` result signals the caller to use the configured default deadband.
    """
    if site_id not in meta.index:
        raise KeyError(f"Site {site_id!r} not found in metadata.")

    if "sigma_t" not in meta.columns:
        return None
    sigma = meta.loc[site_id, "sigma_t"]
    if sigma is None or (isinstance(sigma, float) and math.isnan(sigma)):
        return None
    return float(sigma)
