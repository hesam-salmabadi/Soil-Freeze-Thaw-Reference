# CLAUDE.md

Guidance for Claude Code (and humans) working in this repository.

## Project

**SoFTeR** (Soil Freeze–Thaw Reference) turns heterogeneous in-situ soil
temperature + moisture-probe records into a harmonized, non-binary soil
freeze–thaw reference dataset. It automates the Soil Freezing Characteristic
Curve (SFCC) framework of Salmabadi et al. (2026, *The Cryosphere*,
doi:10.5194/tc-20-1635-2026): instead of a binary 0 °C threshold, each
site–cycle is characterized by a continuous freezing probability
(unfrozen / transitional / frozen) with propagated uncertainty.

The SFCC is fit in **sensor-output-vs-temperature space** — whatever quantity a
probe reports (permittivity, frequency, count, travel/output time). The project
is a pre-release scaffold; it is being built one pipeline stage at a time.

## Commands

Python ≥ 3.10. Environment via conda/mamba or pip.

```bash
# Environment
mamba env create -f environment.yml && mamba activate softer   # or use pip
pip install -e .

# Tests (preferred)
pytest -q

# If an editable install isn't available, run tests against the src layout:
PYTHONPATH=src pytest -q

# Byte-compile check
python -m compileall src/softer

# Lint / format (dev extras)
ruff check src tests
black src tests
```

The CLI entry point is `softer` (`softer.cli:cli`), currently a stub.

## Architecture

`src/` layout, package `softer`. The pipeline mirrors the paper's Materials &
Methods (§2). Modules that are **implemented** are marked; the rest are
docstring-only scaffolds awaiting their stage.

```
src/softer/
  config.py            YAML config loader (typed dataclasses)          [implemented]
  metadata.py          per-site metadata CSV loader (sigma_t)          [implemented]
  sensors.py           sensor -> raw_output registry                   [implemented]
  cli.py, pipeline.py  entry point + stage orchestration               [scaffold]
  io/
    schema.py          common tidy schema                             [scaffold]
    adapters/          per-sensor raw -> schema (HydraProbe/TEROS12/CS616) [scaffold]
  preprocess/          §2.2
    cycles.py          freeze/thaw cycle labeling                      [implemented]
    coverage.py        per-cycle usability / fit-readiness gate        [implemented]
    permittivity.py    raw output -> effective permittivity (eps_eff)  [scaffold]
    binning.py         0.1 degC temperature binning                    [scaffold]
    quality_control.py, resample.py                                    [scaffold]
  model/               §2.3  sfcc, fit, uncertainty, filtering, pooling [scaffold]
  postprocess/         §2.4  probability, classification                [scaffold]
  ancillary/           §2.1.2 era5, ims_snow, soil_grids               [scaffold]
  package/             CF-NetCDF writer + station metadata             [scaffold]
configs/
  config.example.yaml          cross-site constants (see below)
  site_metadata.example.csv    per-site facts: site_id, sensor, sigma_t
tests/                         pytest suite for the implemented modules
```

### Configuration split (important)
- **`configs/*.yaml`** holds values **constant across sites** (windows, durations,
  thresholds, the usability band). Loaded by `softer.config.load_config`.
- **`configs/site_metadata*.csv`** holds **per-site** facts — above all `sigma_t`,
  the soil-temperature uncertainty. `sigma_t` is the **single source of truth** for
  the freeze/thaw deadband; when a site leaves it blank, code falls back to
  `cycle_detection.default_deadband`. Do not reintroduce a sensor-based sigma lookup.

## What the implemented stage does (freeze/thaw labeling)

`softer.preprocess.cycles.label_freeze_thaw_cycles` takes a QA'd per-site series
and adds label columns; only timestamp + soil temperature drive the logic, all
other columns pass through untouched. **No interpolation** — detection runs on the
samples present, so it is insensitive to missing data.

- **Deadband** `h = 2 * sigma_t` (the paper's ±2σ_T). A freezing cycle exists only
  when the soil truly drops below −h; near-zero and small diurnal wiggles are absorbed.
- **Freezing years** run Aug 1 → next Aug 1 (configurable). Detection is independent
  per year; `cycle_id` resets per year.
- **Normal winters** (soil reaches warm before its first cold) use multi-cycle
  detection, keeping every freezing cycle. **Cold-start winters** (begin frozen) are
  flagged and split at the coldest point. Mild / permafrost years set `year_status`.
- **`min_duration`** (default 2 days) absorbs too-short legs into their neighbour.
- **`max_gap`** (default 2 days) stops a cycle spanning a data outage.

Added columns: `freezing_year`, `cycle_id`, `cycle_phase` (1 freezing / 0 thawing /
−1 none), `year_status` (0 normal / 1 never-freezing / 2 always-frozen).

`softer.preprocess.coverage.assess_usability` then adds a boolean `usable` column
and a per-cycle report (writable as a sidecar CSV). Usability is decided on
**temperature-axis coverage** of a freezing cycle: no hole larger than
`max_temp_gap` across the transition band, and it must reach `min_reach`. Time
coverage is a diagnostic flag, not a hard gate.

## Domain glossary
- `sigma_t` — instrument soil-temperature uncertainty; sets the deadband `h = 2σ_T`.
- `EDC` / `eps_eff` — effective (bulk) permittivity; `RDC`/`IDC` are real/imaginary
  parts; `sqrtEDC` = √EDC (the paper fixes the mixing exponent α = 0.5).
- `T_f`, `b` — SFCC freezing onset and transition-sharpness parameters (later stage).
- `eps_int` / `eps_res` — pre-freezing (initial) and residual permittivity.
- freezing / thawing cycle — a cooling leg (peak→trough) / warming leg (trough→peak).

## Conventions
- Keep compute pure and side-effect-free; **plotting stays out of the compute path**
  (optional `plot_*` helpers only). No file I/O inside core functions except the
  explicit `*_csv` / `write_*` wrappers.
- Prefer configurable parameters over hard-coded thresholds; surface new knobs in
  `configs/config.example.yaml` and `softer.config`.
- Add/extend tests in `tests/` for any behavior change; run `pytest` before committing.
- This code is developed collaboratively — do not add comments attributing authorship
  to any tool. Write it as ordinary project code.

## Status / roadmap
Implemented: config + metadata loaders, freeze/thaw cycle labeling, usability gate.
Next stages (scaffolded): eps_eff conversion (`preprocess/permittivity.py`,
`calculate_epsilon_eff`, Mironov), 0.1 °C binning (`preprocess/binning.py`), then
SFCC fit / bootstrap / Monte Carlo / PyMC pooling (`model/`), postprocessing, and
CF-NetCDF packaging.
