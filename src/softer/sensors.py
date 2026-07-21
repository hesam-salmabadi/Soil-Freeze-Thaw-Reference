"""Sensor registry: measurement principles, raw outputs, and uncertainties.

Central place for the per-sensor facts used throughout the pipeline
(Salmabadi et al., 2026, Sect. 2.1.1 and Appendix A):

    - HydraProbe (Stevens Water)  — measures soil temperature + moisture.
    - TEROS12    (METER Group)    — measures soil temperature + moisture.
    - CS616      (Campbell Sci.)  — outputs an oscillation period/count; needs a
                                    separate CS109SS-L temperature sensor.

Each sensor entry captures:
    - ``raw_output`` : the native quantity the probe reports (permittivity,
      frequency, period/count, ...), preserved through ingest.
    - ``sigma_t``    : the instrument-specific soil-temperature uncertainty
      [degC]. This drives the +/-2*sigma_t deadband used for freeze/thaw cycle
      detection, the Monte Carlo temperature perturbation, and freezing-
      probability propagation.

NOTE: the ``sigma_t`` values below are provisional placeholders. They MUST be
reconciled with the sensor-uncertainty table in Appendix A of the paper before
these numbers are used for a published run. Callers can always override the
deadband directly (see ``softer.preprocess.cycles.label_freeze_thaw_cycles``).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SensorSpec:
    """Static facts about a soil-moisture probe used in the study."""

    name: str
    raw_output: str
    sigma_t: float  # instrument-specific soil-temperature uncertainty [degC]


# Provisional sigma_t values — confirm against Appendix A before publication.
_REGISTRY: dict[str, SensorSpec] = {
    "HydraProbe": SensorSpec("HydraProbe", "real_dielectric_permittivity", 0.5),
    "TEROS12": SensorSpec("TEROS12", "raw_counts_permittivity", 0.5),
    "CS616": SensorSpec("CS616", "oscillation_period_us", 0.2),
}


def get_sensor(name: str) -> SensorSpec:
    """Return the :class:`SensorSpec` for ``name`` (case-insensitive-ish lookup)."""
    try:
        return _REGISTRY[name]
    except KeyError:
        # Fall back to a case-insensitive match before giving up.
        for key, spec in _REGISTRY.items():
            if key.lower() == str(name).lower():
                return spec
        known = ", ".join(sorted(_REGISTRY))
        raise KeyError(f"Unknown sensor {name!r}. Known sensors: {known}.") from None


def get_sigma_t(name: str) -> float:
    """Return the instrument-specific soil-temperature uncertainty sigma_t [degC]."""
    return get_sensor(name).sigma_t
