"""Freeze/thaw cycle detection and water-intrusion screening (Sect. 2.2).

Cycle detection based on soil-temperature and eps_eff trends:
    - Freezing cycle: soil temperature decreasing until it reaches its minimum.
    - Thawing cycle : from that minimum until temperature rises above 0 degC.
    - Fluctuations within +/-2*sigma_T of 0 degC (sigma_T = instrument-specific
      temperature uncertainty) are ignored to avoid spurious transient cycles.
    - "Never frozen": if during a freezing cycle T never drops below -sigma_T and
      eps_eff is essentially unchanged, the site/cycle is flagged never-frozen
      (retained for monitoring, but not curve-fit).

Constant-water-content assumption / water-intrusion screening:
    - Total water content is assumed constant through a cycle. Sudden surges in
      eps_eff indicate water entering the system (e.g. snowmelt) and violate the
      assumption; such cycles are excluded.
    - This holds for freezing cycles but generally fails for thawing cycles, so
      the study focuses exclusively on freezing cycles (thawing excluded).

TODO: implement ``detect_cycles(df, sigma_t)``, ``flag_never_frozen(...)``,
and ``screen_water_intrusion(...)``; return freezing cycles only.
"""
