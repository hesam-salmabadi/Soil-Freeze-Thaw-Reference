"""Run configuration loading and validation.

Parses the top-level YAML config (e.g. ``configs/config.example.yaml``) into typed
config objects consumed by the preprocessing stages. Values here are constant
across sites; per-site facts (sigma_t) live in the site-metadata CSV
(:mod:`softer.metadata`).
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _month_day(value) -> tuple[int, int]:
    """Parse ``"MM-DD"`` (or an already-parsed ``(m, d)``) into ``(month, day)``."""
    if isinstance(value, (tuple, list)):
        return int(value[0]), int(value[1])
    month, day = str(value).split("-")
    return int(month), int(day)


def _reject_unknown(cls_name: str, data: dict, allowed: set[str]) -> None:
    extra = set(data) - allowed
    if extra:
        raise KeyError(f"Unknown config key(s) {sorted(extra)} for {cls_name}.")


@dataclass(frozen=True)
class CycleDetectionConfig:
    year_start: tuple[int, int] = (8, 1)
    min_duration: str = "2D"
    max_gap: str = "2D"

    @classmethod
    def from_dict(cls, data: dict) -> "CycleDetectionConfig":
        _reject_unknown(cls.__name__, data, {"year_start", "min_duration", "max_gap"})
        base = cls()
        return cls(
            year_start=_month_day(data.get("year_start", base.year_start)),
            min_duration=data.get("min_duration", base.min_duration),
            max_gap=data.get("max_gap", base.max_gap),
        )


@dataclass(frozen=True)
class UsabilityConfig:
    warm_edge: float = 2.0
    cold_edge: float = -2.0
    min_reach: float = -1.0
    max_temp_gap: float = 0.5
    value_col: str | None = None

    @classmethod
    def from_dict(cls, data: dict) -> "UsabilityConfig":
        allowed = {"warm_edge", "cold_edge", "min_reach", "max_temp_gap", "value_col"}
        _reject_unknown(cls.__name__, data, allowed)
        base = cls()
        return cls(
            warm_edge=float(data.get("warm_edge", base.warm_edge)),
            cold_edge=float(data.get("cold_edge", base.cold_edge)),
            min_reach=float(data.get("min_reach", base.min_reach)),
            max_temp_gap=float(data.get("max_temp_gap", base.max_temp_gap)),
            value_col=data.get("value_col", base.value_col),
        )


@dataclass(frozen=True)
class CriticalWindow:
    start: tuple[int, int] = (10, 1)
    end: tuple[int, int] = (3, 1)

    @classmethod
    def from_dict(cls, data: dict) -> "CriticalWindow":
        _reject_unknown(cls.__name__, data, {"start", "end"})
        base = cls()
        return cls(
            start=_month_day(data.get("start", base.start)),
            end=_month_day(data.get("end", base.end)),
        )


@dataclass(frozen=True)
class CoverageConfig:
    critical_window: CriticalWindow = field(default_factory=CriticalWindow)
    min_coverage: float = 0.80
    max_gap: str = "21D"
    action: str = "flag"

    @classmethod
    def from_dict(cls, data: dict) -> "CoverageConfig":
        allowed = {"critical_window", "min_coverage", "max_gap", "action"}
        _reject_unknown(cls.__name__, data, allowed)
        base = cls()
        return cls(
            critical_window=CriticalWindow.from_dict(data.get("critical_window", {})),
            min_coverage=float(data.get("min_coverage", base.min_coverage)),
            max_gap=data.get("max_gap", base.max_gap),
            action=data.get("action", base.action),
        )


@dataclass(frozen=True)
class Config:
    cycle_detection: CycleDetectionConfig = field(default_factory=CycleDetectionConfig)
    usability: UsabilityConfig = field(default_factory=UsabilityConfig)
    coverage: CoverageConfig = field(default_factory=CoverageConfig)


def load_config(path: str) -> Config:
    """Load and validate a YAML config file into a :class:`Config`."""
    import yaml

    with open(path) as handle:
        raw = yaml.safe_load(handle) or {}

    _reject_unknown("config", raw, {"cycle_detection", "usability", "coverage"})
    return Config(
        cycle_detection=CycleDetectionConfig.from_dict(raw.get("cycle_detection", {})),
        usability=UsabilityConfig.from_dict(raw.get("usability", {})),
        coverage=CoverageConfig.from_dict(raw.get("coverage", {})),
    )
