"""Non-binary soil-state classification (Sect. 2.4).

Maps the continuous freezing probability (:mod:`softer.postprocess.probability`)
onto three physically meaningful states, replacing the binary 0 degC threshold:

    - unfrozen
    - transitional (liquid water and ice coexist)
    - frozen

TODO: implement ``classify(p_frozen)`` with the probability thresholds used to
delimit the three states, returning per-timestamp categorical labels.
"""
