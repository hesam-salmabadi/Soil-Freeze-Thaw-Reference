"""Run configuration loading and validation.

Parses the top-level YAML config (e.g. ``configs/config.example.yaml``) and the
per-network YAML files under ``configs/networks/``, and exposes typed config
objects consumed by :mod:`softer.pipeline`.

A network config declares which adapter to use, the sensor type, file paths /
glob patterns for raw records, and any network-specific overrides. Adding a new
data source should require only a new adapter plus a new network YAML — nothing
else in the pipeline changes.

TODO: define config dataclasses/schema and a ``load_config`` function.
"""
