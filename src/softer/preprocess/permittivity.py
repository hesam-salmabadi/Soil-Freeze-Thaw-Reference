"""Raw sensor output -> effective permittivity (eps_eff), sensor-specific (Sect. 2.2).

Each sensor reports a different native quantity, so the conversion to bulk
effective permittivity is sensor-dependent (details in Appendix A):

    - HydraProbe : manufacturer/real-permittivity output.
    - TEROS12    : raw counts -> permittivity via the manufacturer calibration.
    - CS616      : oscillation period/count -> permittivity.

Absolute eps_eff differences between probes (e.g. from differing operating
frequencies) do not affect soil-state monitoring, because the SFCC is fit
independently per site and cycle (Sect. 2.3.2).

TODO: implement per-sensor conversion functions dispatched by sensor type,
using coefficients from :mod:`softer.sensors`.
"""
