"""Bootstrap and Monte Carlo uncertainty propagation (Sect. 2.3.2).

Block bootstrap:
    - Resample in situ (T_soil, eps_eff) pairs 1000 times.
    - Resample within temperature blocks (with replacement) so the natural
      distribution of eps_eff across the temperature range is preserved while
      introducing between-iteration variability.
    - Refit the SFCC each iteration to obtain parameter means and standard
      deviations for {eps_int, eps_res, b, T_f}.

Monte Carlo (N = 15000 per hourly observation):
    - T_sim ~ Normal(T_obs, sigma_T^2) with sensor-specific sigma_T.
    - Each parameter M in {T_f, eps_int, eps_res, b} sampled independently from
      Normal(mu_M, sigma_M^2) with (mu_M, sigma_M) from the bootstrap.
    - Discard draws violating eps_int < eps_res (physical constraint).
    - Ensemble-mean the modeled eps_fitted per point.

TODO: implement ``block_bootstrap(...)`` and ``monte_carlo(...)``.
"""
