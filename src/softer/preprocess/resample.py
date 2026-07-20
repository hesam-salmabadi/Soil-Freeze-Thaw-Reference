"""Temporal harmonization: hourly resampling and gap interpolation (Sect. 2.2).

Standardizes every record onto a uniform hourly grid with a continuous time
index, then fills gaps by linear interpolation. Interpolated values are flagged
and used *only* for categorical soil-state labeling — they must be excluded from
all curve-fitting analyses to preserve the integrity of model-derived
parameters.

TODO: implement ``resample_hourly(df)`` and ``interpolate_gaps(df)`` with
provenance flagging so fitting stages can drop interpolated points.
"""
