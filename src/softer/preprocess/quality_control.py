"""Quality control of in situ records (Sect. 2.1.1 / 2.2).

Removes physically implausible values before any downstream processing:
    - soil temperature outside [-60, 60] degC,
    - effective permittivity outside [1, 90].

TODO: implement ``apply_qc(df)`` returning the cleaned frame plus quality flags.
"""
