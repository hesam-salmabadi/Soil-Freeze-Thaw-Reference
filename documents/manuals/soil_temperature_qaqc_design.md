# QA/QC Pipeline Design — Multi-Depth Soil Temperature

**Scope.** Automated, sequential quality control for near-surface soil temperature time series. Target depths: 0 cm (surface/skin, optional), **5 cm (primary, near-always present)**, 10 cm (optional). Native resolution varies by site: 30-minute, hourly, or 3-hourly.

**Design basis.** ISMN (Dorigo et al. 2013, geophysical + spectrum/shape categories), IOOS QARTOD (in-situ temperature test set and ordinal flag scheme), NEON (Taylor & Loescher 2013 plausibility tests, sigma/delta persistence, finalQF + quality metrics). Adaptations for a shallow (≤10 cm) stack and for near-0 °C phase-change plateaus are called out inline.

**Status.** Buildable specification. Numeric values are **starting defaults**; per-site, per-depth thresholds should be re-derived from each station's own climatology and sensor noise floor before production.

---

## 1. Conventions

1. **Thresholds are defined in physical time, not sample count.** All window lengths and rate limits are specified in hours/°C-per-hour and converted to sample counts per site using the native interval `Δt`. This keeps one config valid across 30-min / hourly / 3-hourly sites.
2. **Raw is immutable.** Tests write to companion flag columns only; the measured value is never modified or deleted.
3. **Per-test flags are retained separately** (one column per test) and only rolled up at the final aggregation stage. Provenance must survive so that a physically-valid plateau flagged by persistence can be reviewed/overridden rather than silently discarded.
4. **Shape/temporal tests run on contiguous segments.** Missing/structural failures are masked first; spike, step, and persistence tests do not span gaps and emit `NOT_EVALUATED` at segment boundaries.

### Depth bands

| Band | Nominal | Character | Present |
|------|---------|-----------|---------|
| `SURF` | 0 cm | Skin/surface; atmospherically coupled; largest & fastest excursions | Some sites |
| `SHAL` | 5 cm | Primary channel; damped but responsive | Near-always |
| `SUB`  | 10 cm | Most conductive/smoothest of the three | Some sites |

Rationale for treating 0/5/10 cm as one shallow regime: the span (~10 cm) is at or below the diurnal damping depth of typical moist mineral soil (~10–15 cm), so the three channels are strongly coupled rather than spanning distinct thermal regimes. This drives the profile-stage adaptations in §7.

---

## 2. Flag schema (QARTOD-ordinal, non-destructive)

Each test emits one ordinal flag per timestamp:

| Value | Meaning | Downstream action |
|------:|---------|-------------------|
| `1` | PASS | Use |
| `2` | NOT_EVALUATED | Test skipped/masked or could not run (gap boundary, freeze-thaw mask, sensor absent) — **distinct from PASS** |
| `3` | SUSPECT | Retain, downgrade; eligible for human review |
| `4` | FAIL | Reject for science use |
| `9` | MISSING | No datum |

- **Per-test columns:** `qf_<test>_<depth>` (e.g. `qf_spike_5cm`).
- **Summary flag** `qf_summary_<depth>` = worst-of (max) across that depth's per-test flags, with `9` dominating when the value is missing.
- **Publication:** emit as CF-conventions ancillary variables (`flag_values` / `flag_meanings`) alongside each temperature variable. Document every threshold in per-site metadata (ESSD reproducibility expectation).

---

## 3. Processing order

```
T0  Structural        (syntax, missing, gap)
T1  Gross/sensor range (static physical bounds)
T2  Climatology range  (seasonal, per-depth, per-DOY)
T3  Step / rate-of-change (temporal derivative)
T4  Spike              (robust local outlier)
T5  Persistence        (T5a exact-repeat  +  T5b low-variance, freeze-thaw gated)
T6  Profile consistency (opportunistic, ≥2 depths, forcing-gated)
T7  Aggregation        (summary flag + quality metrics)
```

Cheap/deterministic first; shape-based on contiguous segments; cross-depth last. Each stage is depth-aware; the freeze-thaw logic in §6 is cross-cutting and consumed by T4/T5.

---

## 4. Stage specifications (T0–T5)

### T0 — Structural

Purpose: syntax validity, missing-value detection, gap accounting.

```
if value is null/NaN/sentinel:      qf_struct = 9 (MISSING)
elif not parseable/finite:          qf_struct = 4 (FAIL)
elif gap_length ≥ GAP_FAIL_HOURS:   mark surrounding run NOT_EVALUATED for T3–T6
else:                               qf_struct = 1
```
Params: `GAP_FAIL_HOURS = 6 h` (site-tunable). Output: segment boundaries used by T3–T6.

---

### T1 — Gross / sensor range (static)

Purpose: catch electrical faults and physically impossible values.

```
if value < GROSS_MIN or value > GROSS_MAX:   qf_range = 4 (FAIL)
elif value < SENSOR_MIN or > SENSOR_MAX:     qf_range = 4 (FAIL)
else:                                        qf_range = 1
```

| Param | SURF | SHAL | SUB |
|-------|-----:|-----:|----:|
| `GROSS_MIN` (°C) | −45 | −40 | −40 |
| `GROSS_MAX` (°C) | +70 | +55 | +50 |

Surface bounds are **wider**: a skin sensor can exceed air temperature under insolation and undershoot on clear nights. `SENSOR_*` from manufacturer spec.

---

### T2 — Climatology range (seasonal, dynamic)

Purpose: values inside absolute range but implausible for depth + time of year.

```
env = climatology[depth][DOY_window]            # multi-year, ±15-day DOY window
if value < env.p0_5 or value > env.p99_5:       qf_clim = 3 (SUSPECT)
elif |value − env.mean| > K_CLIM · env.sd:      qf_clim = 3
else:                                           qf_clim = 1
if climatology unavailable (short record):      qf_clim = 2 (NOT_EVALUATED)
```
Params: DOY window `±15 d`; `K_CLIM = 4`; envelope from empirical 0.5/99.5 percentiles.
Depth note: build the envelope from **each depth's own** climatology — the seasonal amplitude narrows and lags with depth (SUB narrower than SURF). Do not share bounds across depths.

---

### T3 — Step / rate-of-change

Purpose: unphysically fast change between consecutive samples.

```
roc = |V_i − V_{i-1}| / Δt_hours            # °C per hour
if roc > ROC_FAIL[depth]:      qf_step = 4
elif roc > ROC_SUSPECT[depth]: qf_step = 3
else:                          qf_step = 1
```

| °C/hr | SURF | SHAL | SUB |
|-------|-----:|-----:|----:|
| `ROC_SUSPECT` | 8 | 4 | 3 |
| `ROC_FAIL`    | 12 | 6 | 5 |

Expressed per hour so 30-min / hourly / 3-hourly sites share one config. A jump that is a spike over 30 min can be physical over 3 h — the per-hour normalization handles this. Surface limits are loosest (rain wetting, cloud/insolation transitions are real).

---

### T4 — Spike (robust, Hampel/MAD)

Purpose: local outliers. Preferred over QARTOD neighbor-difference because MAD is robust to the masking problem and to bad neighbors.

```
W  = centered window of SPIKE_WIN_HOURS
m  = median(W)
s  = 1.4826 · MAD(W)
s  = max(s, SIGMA_NOISE[depth])              # NOISE FLOOR — critical, see §6
dev = |V_i − m|
if dev > N_SPIKE_FAIL[depth]    · s: qf_spike = 4
elif dev > N_SPIKE_SUSPECT[depth]· s: qf_spike = 3
else:                                 qf_spike = 1
```

| Param | SURF | SHAL | SUB |
|-------|-----:|-----:|----:|
| `SPIKE_WIN_HOURS` | 3 | 3 | 3 |
| `N_SPIKE_SUSPECT` | 4 | 3 | 3 |
| `N_SPIKE_FAIL`    | 6 | 5 | 4 |
| `SIGMA_NOISE` (°C) | 0.10 | 0.08 | 0.05 |

Depth note: soil is a thermal low-pass filter, so genuine geophysical spikes are effectively absent at 5/10 cm — keep those tight. This does **not** apply to a surface sensor, whose rapid changes are often real — hence looser `N_*` for SURF.
Freeze-thaw note: on an isothermal plateau MAD→0; without the `SIGMA_NOISE` floor, ordinary noise reads as spikes and shreds the zero-curtain. The floor is mandatory.

---

### T5 — Persistence (split test)

The classic sigma/delta/flat-line persistence test flags near-constant runs as "stuck sensor." At these depths that behavior collides directly with the physical near-0 °C plateau. Persistence is therefore split into two sub-tests with different activation.

#### T5a — Exact-repeat stuck test (always on, value-independent)

```
run = count of consecutive bit-identical (or same-quantized) values
if run_hours ≥ REPEAT_FAIL_HOURS:      qf_persist_repeat = 4
elif run_hours ≥ REPEAT_SUSPECT_HOURS: qf_persist_repeat = 3
else:                                  qf_persist_repeat = 1
```
Params: `REPEAT_SUSPECT_HOURS = 6`, `REPEAT_FAIL_HOURS = 12`.
Rationale: a genuine plateau still carries small real fluctuation + digitization noise; a stuck logger repeats an **identical** value. This detects true faults even at 0 °C without touching the curtain, and carries proportionally more of the load here because T5b is masked so often (§6).

#### T5b — Low-variance test (freeze-thaw gated)

```
if |V_i − T_FREEZE| < DELTA_PC[soil]:        qf_persist_var = 2 (NOT_EVALUATED)   # curtain regime
else:
    s_win = SD over PERSIST_WIN_HOURS
    if s_win < VAR_MIN_FAIL:      qf_persist_var = 4
    elif s_win < VAR_MIN_SUSPECT: qf_persist_var = 3
    else:                         qf_persist_var = 1
```

| Param | Default | Note |
|-------|--------:|------|
| `T_FREEZE` (°C) | 0.0 | freezing point |
| `DELTA_PC` (°C) | 0.5 | → 1.0 for saline / fine-textured soils (freezing-point depression) |
| `PERSIST_WIN_HOURS` | 24 | must exceed longest plausible physical plateau at this depth |
| `VAR_MIN_SUSPECT` (°C) | 0.05 | near noise floor |
| `VAR_MIN_FAIL` (°C) | 0.02 | |

Shallow adaptation: near-surface zero-curtains are **short** (hours to a few days at 5 cm, not weeks), and diurnal freeze-thaw *cycling* across 0 °C is common at 5 cm / pervasive at 0 cm in shoulder seasons. The value-gate therefore fires frequently and repeatedly — this is correct; those crossings are the Frozenness-Index signal, not stuck sensors. Set `PERSIST_WIN_HOURS` long enough to clear a genuine curtain but rely on T5a + T6 corroboration inside the phase-change window.

---

## 5. Combined persistence flag

```
qf_persist = max(qf_persist_repeat, qf_persist_var)
# but a T6 profile-corroborated plateau downgrades a T5b SUSPECT to NOT_EVALUATED (see §6/§7)
```

---

## 6. Cross-cutting: freeze-thaw / zero-curtain handling

Physics: during freeze-back or thaw, latent heat of fusion (~334 kJ kg⁻¹) is exchanged at the phase front, pinning temperature near the freezing point (slightly below 0 °C at the freezing-point depression) for extended, genuinely-flat periods. This is the retrieval signal, not a fault. Three mechanisms protect it:

1. **Value-gate (T5b):** suppress the variance-based stuck test inside `|T − T_FREEZE| < DELTA_PC`.
2. **Noise floor (T4):** floor the spike scale estimate at `SIGMA_NOISE` so a flat plateau does not manufacture spikes.
3. **Profile corroboration (T6):** a real curtain propagates **sequentially** down the profile (SURF enters/exits before SHAL before SUB), over hours at these depths. A stuck sensor is flat at one depth while neighbors keep varying. Low variance at depth *k* is only treated as suspect if adjacent depths are **not** showing coherent, time-lagged phase-change behavior.

Decision logic feeding the aggregator:

```
if in_freeze_thaw_window(V):
    if profile_shows_sequential_plateau(depth, t):   # T6
        downgrade qf_persist SUSPECT → 2 (NOT_EVALUATED)   # physical, do not reject
    else:
        keep qf_persist_repeat (T5a) authoritative         # stuck fault still catchable
```

---

## 7. T6 — Profile consistency (opportunistic, forcing-gated)

Runs only when ≥2 depths are present. At 5 cm-only sites this stage is absent and weight shifts to T2 (per-depth climatology) and the T5/§6 freeze-thaw logic, which must be self-sufficient.

Shallow-stack expectations (0/5/10 cm): diurnal amplitude ratio ≈ **1 : 0.65 : 0.45** (SURF:SHAL:SUB); phase lag ≈ **1.5 h per 5 cm** (≈1.6 h at 5 cm, ≈3 h at 10 cm vs. surface).

### T6a — Amplitude damping (swap detector)

```
gate: run only in STRONG_FORCING windows
      (surface diurnal amplitude > AMP_GATE, snow-free, unfrozen, clear-sky)
A_d = rolling diurnal amplitude (SD or MAX−MIN) per depth
expect A_SURF ≥ A_SHAL ≥ A_SUB (monotonic damping)
if ordering violated by > AMP_TOL:   qf_profile = 3 (likely mislabeled/swapped depth)
else:                                qf_profile = 1
```
Params: `AMP_GATE = 4 °C`, `AMP_TOL = 15 %`. Gate is required — when everything is damped (winter, snow-insulated, frozen) all amplitudes collapse to the noise floor and the check is uninformative.

### T6b — Phase lag ordering (resolution-limited)

```
if Δt > 1 h:  qf_phase = 2 (NOT_EVALUATED)         # sub-sample lag at 3-hourly; drop test
else:
    lag_k = argmax cross-correlation(depth_k, depth_{k-1}) over rolling window
    expect 0 < lag_SHAL < lag_SUB (monotonic downward)
    if ordering zero/inverted:  qf_phase = 3        # swap / depth-order error
    else:                       qf_phase = 1
```
At ≤hourly resolution this is the **more reliable** swap detector than T6a when forcing is weak, because ordering survives even when amplitude contrast is tiny.

### T6c — Vertical gradient bounds (pair-specific)

```
grad = |T_upper − T_lower| / Δz_cm            # °C per cm
if grad > GRAD_FAIL[pair]:      qf_grad = 4
elif grad > GRAD_SUSPECT[pair]: qf_grad = 3
else:                           qf_grad = 1
```

| °C/cm | 0–5 cm pair | 5–10 cm pair |
|-------|-----------:|-------------:|
| `GRAD_SUSPECT` | 1.5 | 1.0 |
| `GRAD_FAIL`    | 3.0 | 2.0 |

0–5 cm loosest (steep surface gradients are legitimate under strong forcing); 5–10 cm tighter (more conductive, smoother).

**Dropped for this stack:** deep-convergence-to-mean-annual check — none of these depths is a stable deep anchor; all track diurnal + synoptic forcing.

`qf_profile = max(qf_amp, qf_phase, qf_grad)` where evaluated.

---

## 8. T7 — Aggregation & quality metrics

```
qf_summary[depth] = max over { qf_range, qf_clim, qf_step, qf_spike,
                               qf_persist, qf_profile }        # 9 dominates if missing
```
Apply §6 downgrades before the max. For averaged products (e.g. 30-min → hourly), also emit **quality metrics**: percent of constituent samples at each flag level (NEON-style α/β), so aggregate consumers see partial contamination rather than a single collapsed flag.

**Human-review tier:** segments carrying SUSPECT (`3`) on persistence or profile within the freeze-thaw window are queued for review rather than auto-rejected — these are the highest-value / highest-ambiguity data.

---

## 9. Output columns (per depth)

| Column | Type | Notes |
|--------|------|-------|
| `soil_temp_<depth>` | float | raw, immutable |
| `qf_range_<depth>` … `qf_profile_<depth>` | int8 | per-test ordinal flags |
| `qf_persist_repeat_<depth>`, `qf_persist_var_<depth>` | int8 | split persistence retained |
| `qf_summary_<depth>` | int8 | roll-up |
| `qm_pct_suspect_<depth>`, `qm_pct_fail_<depth>` | float | for averaged products |

CF metadata per flag variable: `flag_values = 1,2,3,4,9`; `flag_meanings = "pass not_evaluated suspect fail missing"`. Publish the full threshold table per site/depth as machine-readable config.

---

## 10. Per-site configuration checklist

Before production at a new site:
1. Set `Δt`; convert all `*_HOURS` params to sample counts.
2. Derive T2 climatology envelopes per present depth (needs multi-year record; else `NOT_EVALUATED`).
3. Measure `SIGMA_NOISE` per sensor from a known-stable quiescent period.
4. Set `DELTA_PC` from soil texture/salinity (0.5 default; up to 1.0 fine/saline).
5. Enumerate present depths → enable/disable T6 and its sub-tests (`Δt > 1 h` disables T6b).
6. Confirm 0 cm sensor type (skin vs. buried) → apply SURF-band widening only if surface/skin.

---

## Appendix — Default parameter summary

| Param | SURF | SHAL | SUB | Global |
|-------|-----:|-----:|----:|-------:|
| GROSS_MIN / MAX (°C) | −45 / +70 | −40 / +55 | −40 / +50 | — |
| K_CLIM | — | — | — | 4 |
| ROC_SUSPECT / FAIL (°C/hr) | 8 / 12 | 4 / 6 | 3 / 5 | — |
| SPIKE_WIN (h) | 3 | 3 | 3 | — |
| N_SPIKE_SUSPECT / FAIL | 4 / 6 | 3 / 5 | 3 / 4 | — |
| SIGMA_NOISE (°C) | 0.10 | 0.08 | 0.05 | — |
| REPEAT_SUSPECT / FAIL (h) | — | — | — | 6 / 12 |
| PERSIST_WIN (h) | — | — | — | 24 |
| VAR_MIN_SUSPECT / FAIL (°C) | — | — | — | 0.05 / 0.02 |
| T_FREEZE / DELTA_PC (°C) | — | — | — | 0.0 / 0.5 |
| AMP_GATE (°C) / AMP_TOL | — | — | — | 4 / 15 % |
| GAP_FAIL (h) | — | — | — | 6 |

*All values are starting defaults for tuning against site climatology and sensor noise, not fixed constants.*
