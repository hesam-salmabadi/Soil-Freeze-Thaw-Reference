"""Sensor registry: measurement principles, raw outputs, and uncertainties.

Central place for the per-sensor facts used throughout the pipeline
(Salmabadi et al., 2026, Sect. 2.1.1 and Appendix A):

    - HydraProbe (Stevens Water)  — measures soil temperature + moisture.
    - TEROS12    (METER Group)    — measures soil temperature + moisture.
    - CS616      (Campbell Sci.)  — outputs an oscillation period/count; needs a
                                    separate CS109SS-L temperature sensor.

Each sensor entry should capture:
    - the raw quantity it reports (permittivity, frequency, period/count,
      travel/output time, ...),
    - the manufacturer conversion to effective permittivity (used by
      :mod:`softer.preprocess.permittivity`),
    - the instrument-specific temperature uncertainty sigma_T used for the
      +/-2*sigma_T cycle-detection deadband, the Monte Carlo temperature
      perturbation, and freezing-probability propagation.

TODO: define a SensorSpec dataclass and a registry keyed by sensor name.
"""
