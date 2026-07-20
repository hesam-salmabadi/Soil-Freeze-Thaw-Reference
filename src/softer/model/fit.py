"""SFCC curve fitting and parameter estimation (Sect. 2.3.2).

Non-linear least-squares fit of the SFCC (:mod:`softer.model.sfcc`) per site and
per freezing cycle, using ``scipy.optimize.curve_fit`` with the Trust Region
Reflective (TRF) algorithm and bounded parameters.

Fitting setup from the paper:
    - alpha fixed at 0.5 (converged to bounds without improving fit).
    - Fit only observations with T_soil <= 2 degC (active freezing range).
    - eps_int: initialized as mean eps_eff within sigma_T <= T_soil <= 2 degC
      (mostly unfrozen Zone 1); bounds = observed [min, max] eps_eff in range.
    - eps_res: lower bound 1 (lowest probe range), upper bound < eps_int.
    - b      : init 1.0, lower bound 0.1, no upper bound.
    - T_f    : allowed up to +1 degC (sensor bias — dielectric may register
               freezing before the thermistor reads subzero).

TODO: implement ``fit_cycle(temps, eps, sensor)`` -> params, covariance, R^2.
"""
