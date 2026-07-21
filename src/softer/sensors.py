"""Sensor registry: measurement principles and native raw outputs.

Reference facts about the soil-moisture probes used in the study (Salmabadi et
al., 2026, Sect. 2.1.1 and Appendix A):

    - HydraProbe (Stevens Water)  — measures soil temperature + moisture.
    - TEROS12    (METER Group)    — measures soil temperature + moisture.
    - CS616      (Campbell Sci.)  — outputs an oscillation period/count; needs a
                                    separate CS109SS-L temperature sensor.

Each entry records ``raw_output`` — the native quantity the probe reports
(permittivity, frequency, period/count, ...), which the ε_eff conversion step
consumes. The soil-temperature uncertainty ``sigma_t`` that sets the freeze/thaw
deadband is NOT stored here: it is a per-site value and lives solely in the
site-metadata CSV (:mod:`softer.metadata`).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SensorSpec:
    """Static facts about a soil-moisture probe used in the study."""

    name: str
    raw_output: str


_REGISTRY: dict[str, SensorSpec] = {
    "HydraProbe": SensorSpec("HydraProbe", "real_dielectric_permittivity"),
    "TEROS12": SensorSpec("TEROS12", "raw_counts_permittivity"),
    "CS616": SensorSpec("CS616", "oscillation_period_us"),
}


def get_sensor(name: str) -> SensorSpec:
    """Return the :class:`SensorSpec` for ``name`` (case-insensitive fallback)."""
    try:
        return _REGISTRY[name]
    except KeyError:
        for key, spec in _REGISTRY.items():
            if key.lower() == str(name).lower():
                return spec
        known = ", ".join(sorted(_REGISTRY))
        raise KeyError(f"Unknown sensor {name!r}. Known sensors: {known}.") from None
