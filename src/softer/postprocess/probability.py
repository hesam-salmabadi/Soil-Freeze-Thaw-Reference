"""Freezing probability of seasonally frozen ground (Sect. 2.4.1).

Computes P_frozen (the freezing probability / degree of soil freezing) at each
timestamp by propagating uncertainty from the hierarchical posteriors of T_f and
b:
    - Draw paired posterior samples {(T_f^(s), b^(s))} from the PyMC models,
      restoring b^(s) = exp(b_j^(s)) to the original scale.
    - Perturb soil-temperature observations by sensor uncertainty:
      T_sim^(s) ~ Normal(T_obs, sigma_T^2).
    - Evaluate the normalized SFCC (:mod:`softer.model.sfcc`) for each draw and
      average across Monte Carlo samples to get the mean freezing probability
      per timestamp.

TODO: implement ``freezing_probability(t_obs, posterior, sensor)``.
"""
