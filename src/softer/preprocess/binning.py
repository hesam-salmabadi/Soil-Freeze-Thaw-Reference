"""Balanced-dataset temperature binning (Sect. 2.2).

The distribution of observations across temperature ranges is usually uneven,
which biases curve fitting. To balance it, eps_eff values are averaged within
0.1 degC soil-temperature bins (matching the sensors' temperature resolution).
This both compensates for the uneven distribution and reduces noise from diurnal
temperature fluctuations that are absent in controlled lab settings.

TODO: implement ``bin_by_temperature(df, width=0.1)`` returning per-bin mean
eps_eff (and bin membership counts for diagnostics).
"""
