"""Post-fit filtering of SFCC cycles (Sect. 2.3.2).

Multi-step quality filtering to keep only physically valid, reliable fits:
    - Keep only fall/winter cycles (1 Sep - 1 Mar, the "freezing season");
      exclude short transient spring events.
    - Remove cycles with R^2 < 0.6.
    - Remove visually detected anomalies (e.g. irregular water-content changes
      during freezing).
    - Remove unreliable parameter estimates: boundary values for T_f or b, or
      excessively wide confidence intervals.
    - Remove extreme outliers in T_f and b by retaining only the central 95%
      (2.5th-97.5th percentiles).

TODO: implement ``filter_cycles(fits)`` applying the criteria above.
"""
