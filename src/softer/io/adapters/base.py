"""Adapter interface.

Defines the common contract every source adapter implements so the pipeline can
treat all networks uniformly. An adapter is responsible only for *reading* a
network's raw files and emitting rows in the common schema
(:mod:`softer.io.schema`) — it does not do QC, unit harmonization beyond what is
sensor-native, or cycle detection (those live in :mod:`softer.preprocess`).

TODO: define ``BaseAdapter`` with, e.g.::

    class BaseAdapter(ABC):
        sensor: str
        @abstractmethod
        def read(self, source) -> "pandas.DataFrame": ...
"""
