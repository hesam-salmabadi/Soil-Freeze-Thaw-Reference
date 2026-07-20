# SoFTeR

**So**il **F**reeze–**T**haw **Re**ference.

**A non-binary in situ reference dataset of soil freeze–thaw state and soil moisture
for calibration and validation of remote sensing products.**

SoFTeR automates the permittivity–temperature Soil Freezing Characteristic Curve
(SFCC) framework of Salmabadi et al. (2026, *The Cryosphere*,
[10.5194/tc-20-1635-2026](https://doi.org/10.5194/tc-20-1635-2026)) and applies it
to heterogeneous in situ soil temperature and permittivity / soil moisture records
from multiple monitoring networks. It produces a harmonized reference dataset in
which each site–cycle is characterized by continuous freezing probability
(unfrozen / transitional / frozen), rather than a binary 0 °C threshold, together
with propagated uncertainty. The dataset is intended for evaluation of satellite
freeze–thaw and soil moisture products, and is packaged for publication in
*Earth System Science Data* (ESSD).

## Pipeline

Raw source records are processed through a configurable pipeline:

1. **Ingest** — source adapters convert raw sensor outputs to a common schema.
2. **Preprocess** — raw → effective permittivity (sensor-specific); freeze/thaw
   cycle detection; balanced 0.1 °C temperature binning.
3. **Model** — SFCC fit in permittivity–temperature space (SciPy TRF, bounded);
   block bootstrap; Monte Carlo uncertainty propagation.
4. **Pool** — Bayesian hierarchical partial pooling of `b` and `T_f`
   (network / site random effects) in PyMC.
5. **Postprocess** — freezing probability and soil-state classification.
6. **Package** — CF-compliant NetCDF + station metadata for release.
7. **Cal/val** — optional collocation against remote sensing products.

## Quick start

```bash
mamba env create -f environment.yml
mamba activate softer
pip install -e .
softer run --config configs/config.example.yaml
```

## Repository layout

See `docs/` for the methodology note and data dictionary. Adding a new data source
means writing one adapter under `src/softer/io/adapters/` and one YAML file under
`configs/networks/`; nothing else in the pipeline changes.

## Status

Pre-release scaffold. Not yet a citable dataset.

## License

Code: TBD (see `LICENSE`). Data release: CC-BY-4.0 recommended for ESSD.
