# SoFTeR

**So**il **F**reeze–**T**haw **Re**ference.

**A non-binary in situ reference dataset of soil freeze–thaw state
for calibration and validation of remote sensing products.**

SoFTeR generalizes the Soil Freezing Characteristic Curve (SFCC) framework of
Salmabadi et al. (2026, *The Cryosphere*,
[10.5194/tc-20-1635-2026](https://doi.org/10.5194/tc-20-1635-2026)) and applies it
to heterogeneous in situ soil temperature records paired with the raw output of
soil-moisture probes — whatever quantity a given sensor reports (permittivity,
oscillation frequency, count, travel/output time, output period, etc.) — from
multiple monitoring networks. The SFCC is fit in that sensor-output–versus–temperature
space rather than being restricted to permittivity. It produces a harmonized
reference dataset in which each site–cycle is characterized by continuous freezing
probability (unfrozen / transitional / frozen), rather than a binary 0 °C threshold,
together with propagated uncertainty. The dataset is intended for evaluation of
satellite freeze–thaw products, and is packaged for publication in
*Earth System Science Data* (ESSD).

## Pipeline

Raw source records are processed through a configurable pipeline:

1. **Ingest** — source adapters convert raw sensor outputs to a common schema,
   preserving whatever quantity each probe reports (permittivity, frequency,
   count, travel/output time, etc.).
2. **Preprocess** — sensor-specific handling of the raw probe output; freeze/thaw
   cycle detection; balanced 0.1 °C temperature binning.
3. **Model** — SFCC fit in sensor-output–versus–temperature space (SciPy TRF,
   bounded); block bootstrap; Monte Carlo uncertainty propagation.
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
