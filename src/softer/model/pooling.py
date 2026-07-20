"""Bayesian hierarchical partial pooling of b and T_f (Sect. 2.3.2).

Stabilizes per-cycle estimates of the shape factor b and freezing onset T_f by
sharing information across sites and networks within ecozones (and, separately,
at the biome level), reducing uncertainty while preserving genuine spatial
differences. Implemented in PyMC.

Model structure:
    - b   : model log(b) with network- and site-level random intercepts,
            non-centered parameterization, Student-t likelihood whose
            heteroskedastic scale is the bootstrap SE of log(b) plus a residual.
    - T_f : analogous hierarchy on the original degC scale, with observation-level
            SEs and an additional residual variance.
    - Priors: weakly informative Normal on means; Half-Student-t on variance
      components.
    - Inference: NUTS. Confirm convergence (R-hat <= 1.01, high ESS), stability
      (no/rare divergences), energy (E-BFMI >= 0.5), and PSIS-LOO (all k-hat<=0.7).

TODO: implement ``fit_hierarchy(cycles)`` returning posterior samples of b, T_f
per network/site (and biome-level variant).
"""
