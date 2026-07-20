"""Common in-memory schema shared across pipeline stages.

Every adapter (:mod:`softer.io.adapters`) converts a network's raw records into
this single tidy schema so that all downstream stages are source-agnostic.

Minimal per-observation columns (proposed):
    - ``time``          : timezone-aware UTC timestamp (hourly after resampling)
    - ``network``       : network code (FM, LR, BJ, CP, BT, KN, TV, GR)
    - ``site``          : site identifier within the network
    - ``sensor``        : sensor type (HydraProbe / TEROS12 / CS616)
    - ``t_soil``        : soil temperature [degC]
    - ``raw_output``    : the sensor's native moisture output (as reported)
    - ``eps_eff``       : effective permittivity after sensor-specific conversion
    - ``quality_flag``  : QC / interpolation provenance flag

TODO: define the column contract (dtypes, units) and light validators.
"""
