"""End-to-end pipeline orchestration.

Wires together the processing stages in order, mirroring the methodology of
Salmabadi et al. (2026), Sect. 2:

    1. Ingest      -> softer.io.adapters      (Sect. 2.1.1)
    2. Preprocess  -> softer.preprocess       (Sect. 2.2)
    3. Model       -> softer.model            (Sect. 2.3)
    4. Pool        -> softer.model.pooling    (Sect. 2.3.2)
    5. Postprocess -> softer.postprocess      (Sect. 2.4)
    6. Package     -> softer.package          (release)
    7. Cal/val     -> (optional) collocation against RS products

Each stage is a thin, testable function operating on a common in-memory schema
(:mod:`softer.io.schema`). The orchestrator here is deliberately dumb: it reads
config, calls each stage, and hands intermediate results to the next.

TODO: implement ``run(config)`` that executes the stages above.
"""
