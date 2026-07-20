"""CS616 (Campbell Scientific) adapter — Sect. 2.1.1, Appendix A.

Network: Chapleau (CP). The CS616 outputs an oscillation period/count rather
than permittivity or temperature directly, so this adapter must also pair each
probe with its separate CS109SS-L soil-temperature sensor.

Chapleau specifics:
    - 24 CS616 probes across four 200 m x 200 m plots (4 sites, 24 probes);
      multiple probes per plot are averaged to represent plot conditions.
    - The 30 cm needle is angled at 20 deg to integrate the top 10 cm of soil
      (midpoint depth 5 cm), for comparability with the 5 cm sites.

TODO: implement ``CS616Adapter(BaseAdapter)`` including the temperature-sensor
pairing and per-plot probe averaging.
"""
