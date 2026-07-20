"""Preprocessing stage (Sect. 2.2).

Raw sensor output -> effective permittivity, quality control, uniform hourly
resampling with gap interpolation, freeze/thaw cycle detection, water-intrusion
screening, and 0.1 degC temperature binning to produce a balanced dataset for
curve fitting.
"""
