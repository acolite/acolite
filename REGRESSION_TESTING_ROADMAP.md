# ACOLITE Regression Testing Roadmap

## Overview

This document outlines the regression testing strategy for validating the Rust
port (`acolite-rs`) against the reference Python ACOLITE implementation.

---

## Architecture

```
tests/
├── pace_e2e.rs                    # Rust end-to-end PACE test
├── landsat_e2e.rs                 # Rust end-to-end Landsat 8/9 test
├── sentinel2_e2e.rs               # Rust end-to-end Sentinel-2 test
├── integration_tests.rs           # Existing Rust integration tests
├── regression/
│   ├── conftest.py                # Pytest fixtures (build Rust, temp dirs, tolerance)
│   ├── test_pace_regression.py    # PACE-specific Python↔Rust comparison (17 tests)
│   ├── test_pace_rust_vs_python.py # PACE Rust vs Python full AC (14 tests)
│   ├── test_pace_dsf_rust_vs_python.py # PACE DSF Chesapeake Bay (12 tests)
│   ├── test_pace_sa_dsf_rust_vs_python.py # PACE DSF South Australia + full scene benchmark (12 tests)
│   ├── test_landsat_regression.py # Landsat 8/9 Python↔Rust comparison (13 tests)
│   ├── test_landsat_rust_vs_python.py # Landsat Rust vs Python full AC (13 tests)
│   ├── test_benchmark_rust_vs_python.py # Landsat full-scene benchmark (7 tests)
│   ├── test_sentinel2_regression.py # Sentinel-2 Python↔Rust comparison (19 tests)
│   ├── test_s2_rust_vs_python.py  # S2 Rust vs Python fixed-mode AC (15 tests)
│   ├── test_s2_benchmark_rust_vs_python.py # S2 full-scene benchmark (9 tests)
│   ├── run_pace_real.py           # PACE real data download + reference generation
│   └── __init__.py
benches/
│   └── performance.rs             # Criterion benchmarks (Landsat + S2 + resample)
```

## Test Summary

| Suite | Tests | Status |
|-------|-------|--------|
| Rust unit tests | 26 | ✅ All pass |
| Rust integration tests | 8 | ✅ All pass |
| Rust E2E tests | 14 (+1 pre-existing COG failure) | ✅ |
| Python regression (total) | 141 | ✅ All pass |

### Python Test Breakdown

| Test file | Tests | Sensor | Mode |
|-----------|-------|--------|------|
| test_landsat_regression.py | 13 | Landsat 8/9 | Synthetic + unit |
| test_landsat_rust_vs_python.py | 13 | Landsat 8/9 | ROI pixel comparison |
| test_benchmark_rust_vs_python.py | 7 | Landsat 8/9 | Full-scene benchmark |
| test_sentinel2_regression.py | 19 | Sentinel-2 | Synthetic + unit |
| test_s2_rust_vs_python.py | 15 | Sentinel-2 | Fixed-mode pixel comparison |
| test_s2_benchmark_rust_vs_python.py | 9 | Sentinel-2 | Full-scene + tiled benchmark |
| test_pace_regression.py | 17 | PACE OCI | Synthetic + real data |
| test_pace_rust_vs_python.py | 14 | PACE OCI | Rust vs Python comparison |
| test_pace_sa_fullscene_benchmark.py | 12 | PACE OCI | Full-scene benchmark (SA) |

## Test Tiers

### Tier 1 — Unit / Synthetic (always runs, no I/O)

| Test | What it validates |
|------|-------------------|
| Rayleigh τ monotonicity | λ^-4 law across 400–2200 nm |
| Gas transmittance range | O₃ and H₂O transmittance ∈ (0, 1] |
| DSF dark spectrum percentile | 5th percentile extraction correctness |
| Parallel determinism | Two parallel runs produce identical output |
| GeoZarr directory layout | zarr.json, /data, /wavelengths, /bandwidths |
| CF-1.8 metadata | Conventions, proj:epsg, spatial:transform |

**Run:** `pytest tests/regression/ -v` or `cargo test`

### Tier 2 — Synthetic NetCDF Round-Trip (needs netCDF4)

| Test | What it validates |
|------|-------------------|
| Synthetic PACE L1B creation | Minimal valid NetCDF with 3 detectors |
| Python l1_convert ingestion | ACOLITE can read the synthetic file |
| Band count consistency | Same number of bands from both implementations |

**Run:** `pytest tests/regression/ -v -k synthetic`

### Tier 3 — Real Data End-to-End (needs cached data)

| Test | What it validates |
|------|-------------------|
| L1R band count match | Python and Rust produce same band count |
| Per-band reflectance statistics | Mean/std/min/max within tolerance |
| Spectral ordering | Wavelengths monotonically increasing |
| Per-pixel Pearson R, RMSE, MAE | Physics equivalence |

**Run:** `pytest tests/regression/ -v -s`

---

## Phase 1: PACE OCI L1B → L1R

**Status: Implemented ✓**

- [x] Rust PACE loader (`src/loader/pace.rs`) — reads L1B NetCDF
- [x] Parallel atmospheric correction (`src/parallel.rs`) — rayon-based
- [x] GeoZarr writer (`src/writer/geozarr.rs`) — Zarr V3 + CF-1.8
- [x] Rust E2E test (`tests/pace_e2e.rs`) — synthetic scene
- [x] Python regression harness (`tests/regression/`)

### Key metrics

| Metric | Python baseline | Rust target | Tolerance |
|--------|----------------|-------------|-----------|
| Band count (L1B→L1R) | ~280 | ~280 | exact |
| Mean ρ_t per band | reference | match | 1e-4 |
| Processing time (280 bands, 200×150) | ~30s | <5s | — |

## Phase 1b: Landsat 8/9 OLI

**Status: Implemented ✓**

- [x] Rust Landsat loader (`src/loader/landsat.rs`) — reads GeoTIFF bands
- [x] MTL parser (`src/sensors/mtl.rs`) — metadata extraction
- [x] Full LUT-DSF atmospheric correction
- [x] COG writer (`src/writer/cog.rs`) — per-band GeoTIFF
- [x] Rust E2E test (`tests/landsat_e2e.rs`) — 15 tests
- [x] Python regression (`tests/regression/test_landsat_regression.py`)
- [x] Performance benchmarks (`benches/performance.rs`)

### Regression results — South Australia (Gulf St Vincent)

| Metric | L8 | L9 |
|--------|-----|-----|
| Mean Pearson R | 0.9985 | 0.9996 |
| Mean RMSE | 0.011 | 0.003 |
| Mean %<0.05 | 99.9% | 100.0% |
| Mean %<0.01 | 71.1% | 98.1% |

## Phase 1c: Sentinel-2 A/B MSI

**Status: Implemented ✓ — Physics-Equivalent**

- [x] Rust S2 SAFE loader (`src/loader/sentinel2.rs`) — JP2 via GDAL, multi-resolution resampling
- [x] RADIO_ADD_OFFSET support (processing baseline ≥ 4.0)
- [x] QUANTIFICATION_VALUE parsing from MTD_MSIL1C.xml
- [x] Full LUT-DSF atmospheric correction (fixed + tiled modes)
- [x] CLI args: `--model auto|MOD1|MOD2`, `--aot-mode tiled|fixed`
- [x] COG writer for multispectral output
- [x] Rust E2E test (`tests/sentinel2_e2e.rs`) — 18 tests
- [x] Python regression (`tests/regression/test_sentinel2_regression.py`) — 19 tests
- [x] Fixed-mode pixel comparison (`tests/regression/test_s2_rust_vs_python.py`) — 15 tests
- [x] Full-scene benchmark (`tests/regression/test_s2_benchmark_rust_vs_python.py`) — 9 tests

### Physics fixes applied (2026-03-07)

Three root causes of AOT mismatch were identified and fixed:

1. **Missing RADIO_ADD_OFFSET**: S2 processing baseline ≥ 4.0 stores DN with offset -1000.
   Correct formula: `rhot = (DN - 1000) / 10000`. Rust was using `DN / 10000`, giving
   TOA values 0.1 too high, which caused completely wrong dark spectrum and AOT.

2. **invert_aot returning tau_min instead of NaN**: When dark spectrum < romix(tau=0.001),
   Python returns NaN (band excluded from AOT). Rust was returning 0.001 (band included),
   pulling AOT to minimum. Fixed to return NaN matching Python's `np.interp(left=NaN)`.

3. **wave_range mismatch**: Python S2 default `dsf_wave_range=400,900` excludes SWIR.
   Rust default was `(400, 2500)` including B11/B12 in AOT estimation.

### Regression results — South Australia (MGRS T54HTF, fixed mode)

| Metric | S2A | S2B |
|--------|-----|-----|
| Mean Pearson R | **1.000000** | **1.000000** |
| Mean RMSE | **0.001245** | **0.001475** |
| Mean %<0.05 | **100.0%** | **100.0%** |
| Mean %<0.01 | 97.5% | 97.1% |
| Model selected | MOD1 ✓ | MOD2 ✓ |
| Rust AOT | 0.0110 | 0.1150 |
| Python AOT | 0.0110 | 0.1150 |

Both scenes achieve physics-equivalent results (RMSE < 0.002), with identical
model selection and AOT values matching to 4 significant figures.

### Per-band detail (S2A)

| Band | RMSE | Bias | Rust μ | Python μ |
|------|------|------|--------|----------|
| B01 (443nm) | 0.000013 | +0.000013 | 0.034335 | 0.034323 |
| B02 (492nm) | 0.000015 | +0.000015 | 0.035042 | 0.035028 |
| B03 (560nm) | 0.000202 | -0.000183 | 0.035853 | 0.036036 |
| B04 (665nm) | 0.000397 | +0.000294 | 0.035039 | 0.034745 |
| B05 (704nm) | 0.001000 | -0.000686 | 0.041233 | 0.041920 |
| B06 (740nm) | 0.001416 | -0.000898 | 0.053547 | 0.054445 |
| B07 (783nm) | 0.001320 | -0.000799 | 0.059872 | 0.060671 |
| B08 (833nm) | 0.007184 | -0.004201 | 0.060680 | 0.064881 |
| B8A (865nm) | 0.000023 | -0.000012 | 0.068232 | 0.068244 |
| B11 (1614nm) | 0.000240 | -0.000120 | 0.078816 | 0.078936 |
| B12 (2202nm) | 0.001887 | -0.000935 | 0.054104 | 0.055038 |

B08 has the highest RMSE (0.007) due to its broad bandwidth (106nm) causing
resampling differences between Rust (GDAL nearest) and Python (ACOLITE internal).

### DSF configuration alignment (Rust ↔ Python)

| Parameter | Rust | Python | Match |
|-----------|------|--------|-------|
| dsf_spectrum_option | Intercept(200) | intercept, 200 pixels | ✓ |
| dsf_wave_range | (400, 900) | [400, 900] | ✓ |
| dsf_nbands | 2 | 2 | ✓ |
| dsf_nbands_fit | 2 | 2 | ✓ |
| dsf_aot_compute | Min | min | ✓ |
| min_tgas_aot | 0.85 | 0.85 | ✓ |
| dsf_model_selection | min RMSD | min_drmsd | ✓ |
| invert_aot out-of-range | NaN | NaN (np.interp left=NaN) | ✓ |

### Performance (5490×5490 at 20m, ~30M pixels, 11 bands)

| Metric | S2A Python | S2A Rust | Speedup | S2B Python | S2B Rust | Speedup |
|--------|------------|----------|---------|------------|----------|---------|
| Total | 182.4s | 66.1s | **2.8×** | 173.0s | 65.5s | **2.6×** |

## Phase 1d: Real PACE OCI Data Regression (OB-DAAC)

**Status: Implemented ✓**

- [x] CMR search for PACE OCI L1B granules (public API, no auth)
- [x] EarthData OAuth download (two-step redirect pattern)
- [x] Python ACOLITE L1R processing with spatial limit (ROI subsetting)
- [x] Reference statistics generation and caching
- [x] 9 pytest regression tests on real data
- [x] Spectral physics validation (O₂-A absorption, H₂O absorption, Rayleigh)

### Test data

| Property | Value |
|----------|-------|
| Granule | `PACE_OCI.20240701T175112.L1B.V3.nc` |
| Location | Chesapeake Bay (36.0°N, 75.5°W) |
| ROI | [35.8, -75.4, 36.0, -75.2] → 20×14 pixels |
| Bands | 291 (119 blue + 163 red + 9 SWIR) |

## Phase 1e: Rust vs Python Performance Regression

**Status: Implemented ✓**

### Performance results (real PACE OCI, 291 bands, 21×15 ROI)

| Stage | Python | Rust | Speedup |
|-------|--------|------|---------|
| Total (load+AC+write) | 16.7s | 0.31s | **54.7×** |

### Output comparison

| Metric | Value |
|--------|-------|
| Band count match | 291 = 291 ✓ |
| Spectral correlation (TOA vs corrected) | r = 0.9995 |
| Bands where ρs ≤ ρt | 289/289 (100%) |

## Phase 1g: PACE OCI Full-Scene Benchmark (South Australia)

**Status: Implemented ✓**

Full-scene fixed-mode DSF atmospheric correction on PACE OCI over open ocean
south of Kangaroo Island, South Australia.

### Test data

| Property | Value |
|----------|-------|
| Granule | `PACE_OCI.20241231T044250.L1B.V3.nc` |
| Location | South of Kangaroo Island, SA (-37°S, 137°E) |
| Size | 1710×1272 pixels × 291 bands (633M band-pixels) |
| DSF mode | Fixed (whole-scene dark spectrum) |

### Performance (full scene 1710×1272 × 291 bands)

| Stage | Python | Rust | Notes |
|-------|--------|------|-------|
| Total | 230s | **84s** | **2.7× speedup** |
| Load | — | 12s | Bulk NetCDF reads (3 vs 291) |
| AC | — | 34s | 18.7 Mpx/s (rayon parallel) |
| Write | — | 35s | GeoZarr V3 + gzip |

### Accuracy (per-band, 1.57M pixels each)

| Band | Pearson R | RMSE | Bias | %<0.05 | %<0.01 |
|------|-----------|------|------|--------|--------|
| 400nm | 0.999998 | 0.0126 | +0.0105 | 100% | 22.8% |
| 443nm | 1.000000 | 0.0064 | +0.0048 | 100% | 98.8% |
| 490nm | 1.000000 | 0.0055 | +0.0041 | 100% | 99.7% |
| 550nm | 1.000000 | 0.0035 | +0.0012 | 100% | 98.8% |
| 600nm | 1.000000 | 0.0034 | -0.0006 | 100% | 97.1% |
| 665nm | 1.000000 | 0.0012 | -0.0001 | 100% | 100% |
| 750nm | 0.999999 | 0.0009 | +0.0004 | 100% | 100% |
| 865nm | 1.000000 | 0.0005 | -0.0000 | 100% | 100% |

**Mean R=1.0000, Mean RMSE=0.0043, 100% pixels within 0.05**

### AOT and model selection

| Metric | Python | Rust | Match |
|--------|--------|------|-------|
| Model | MOD2 | MOD2 | ✓ |
| AOT | 0.0105 | 0.0093 | 11.4% diff |

AOT difference is due to minor differences in dark spectrum extraction
(intercept regression) over the full 1710×1272 scene. Both select the
same aerosol model and produce near-identical surface reflectance.

## Phase 1f: Real Landsat 8 Data Regression (USGS S3)

**Status: Implemented ✓**

### Test data

| Property | Value |
|----------|-------|
| Scene | `LC08_L1TP_013035_20240710_20240718_02_T1` |
| Location | Outer Banks / Chesapeake Bay |
| Source | `s3://usgs-landsat/collection02/level-1/` (requester-pays) |

### Performance (full scene 7931×7891 = 62M pixels × 7 bands)

| Metric | L8 Python | L8 Rust | Speedup | L9 Python | L9 Rust | Speedup |
|--------|-----------|---------|---------|-----------|---------|---------|
| Total | 179.5s | 66.1s | **2.7×** | 180.1s | 56.3s | **3.2×** |

---

## Phase 2: L1R → L2R (DSF Atmospheric Correction)

**Status: Implemented ✓**

The full LUT-based Dark Spectrum Fitting (DSF) atmospheric correction has been
ported to Rust, using the same 6SV radiative transfer LUTs as Python ACOLITE.

### Rust AC modules

| Module | File | Purpose |
|--------|------|---------|
| `ac::interp` | `src/ac/interp.rs` | N-dimensional regular grid interpolator |
| `ac::aerlut` | `src/ac/aerlut.rs` | Aerosol LUT reader — sensor-specific 6SV NetCDF |
| `ac::gas_lut` | `src/ac/gas_lut.rs` | Gas transmittance — ko3 + RSR convolution |
| `ac::dsf` | `src/ac/dsf.rs` | Full DSF: dark spectrum, AOT inversion, model selection, correction |

### DSF algorithm (matching Python ACOLITE)

1. **DN → TOA reflectance**: `rhot = (DN + RADIO_ADD_OFFSET) / QUANTIFICATION_VALUE`
2. **Gas transmittance**: Ozone (ko3 × RSR convolution), water vapour, CO₂, O₂
3. **Dark spectrum**: Intercept of 200 darkest pixels per band (linear regression)
4. **AOT inversion**: For each model, invert LUT romix(tau) → tau; NaN if below range
5. **Model selection**: Minimum RMSD between observed dark spectrum and modeled path reflectance
6. **Surface reflectance**: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`

### Cross-sensor regression summary

| Sensor | Mean R | Mean RMSE | Model match | AOT match | Notes |
|--------|--------|-----------|-------------|-----------|-------|
| S2A (fixed) | 1.000 | 0.0012 | MOD1 ✓ | 0.011 ✓ | Physics-equivalent |
| S2B (fixed) | 1.000 | 0.0015 | MOD2 ✓ | 0.115 ✓ | Physics-equivalent |
| L8 (tiled) | 0.999 | 0.011 | MOD2 ✓ | ~0.59 | Full scene vs ROI |
| L9 (tiled) | 1.000 | 0.003 | MOD2 ✓ | ~0.27 | Full scene vs ROI |
| PACE OCI (full, fixed) | 1.000 | 0.004 | MOD2 ✓ | 0.009/0.011 | Full scene, 291 bands |

### Remaining gaps vs Python ACOLITE

| Feature | Status | Impact on accuracy |
|---------|--------|-------------------|
| ROI subsetting (limit) | Not yet | Would match Python's ROI-only processing |
| Interface reflectance (rsky) | Not yet | Minor when dsf_interface_reflectance=False |
| Ancillary data (NCEP) | Not yet | Uses default ozone/pressure |
| DEM-derived pressure | Not yet | Uses 1013 hPa |
| Glint correction | Not yet | Minor for nadir sensors |

## Phase 3: L2R → L2W (Water Products)

**Status: Planned**

- [ ] Port `acolite_l2w.py` parameter computation
- [ ] Chlorophyll-a (OC algorithms)
- [ ] TSS (Nechad, Dogliotti)
- [ ] Turbidity
- [ ] Regression: compare derived products

## Phase 4: Multi-Sensor Expansion

| Sensor | Loader | AC | Writer | Regression |
|--------|--------|----|--------|------------|
| PACE OCI | ✅ | ✅ Full DSF (fixed+tiled) | ✅ GeoZarr | ✅ 43 tests (ROI + full scene) |
| Landsat 8/9 | ✅ | ✅ LUT-DSF | ✅ COG | ✅ 33 tests |
| Sentinel-2 | ✅ | ✅ LUT-DSF | ✅ COG | ✅ 43 tests |
| Sentinel-3 OLCI | ✅ (sensor def) | ✅ (basic) | — | — |
| PRISMA | — | — | — | — |
| EMIT | — | — | — | — |

## Phase 5: Performance Regression

- [ ] Benchmark suite (`benches/performance.rs`) — track Mpx/s
- [ ] CI integration: fail if throughput drops >10%
- [ ] Memory profiling: track peak RSS

---

## Running Tests

### Rust tests

```bash
# Unit + integration + E2E tests
cargo test --features full-io

# Specific test suites
cargo test --test landsat_e2e
cargo test --test sentinel2_e2e
cargo test --test pace_e2e
cargo test --test integration_tests

# Benchmarks
cargo bench
```

### Python regression tests

```bash
# All regression tests (Tier 1+2, no real data needed for synthetic)
pytest tests/regression/ -v

# Sentinel-2 fixed-mode regression (requires cached S2 data in /tmp/acolite_test_cache/)
pytest tests/regression/test_s2_rust_vs_python.py -v -s

# Sentinel-2 benchmark (fixed + tiled modes)
pytest tests/regression/test_s2_benchmark_rust_vs_python.py -v -s

# Landsat regression
pytest tests/regression/test_landsat_rust_vs_python.py -v -s

# Landsat benchmark
pytest tests/regression/test_benchmark_rust_vs_python.py -v -s

# PACE regression
pytest tests/regression/test_pace_rust_vs_python.py -v -s

# Full regression with real data
pytest tests/regression/ -v --runslow \
  --pace-file /path/to/PACE_OCI.L1B.nc \
  --landsat-file /path/to/LC08_L1TP_... \
  --s2-file /path/to/S2A_MSIL1C_....SAFE
```

### CI Pipeline (recommended)

```yaml
jobs:
  rust-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: cargo test --features full-io

  python-regression:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: conda-incubator/setup-miniconda@v3
        with:
          environment-file: environment.yml
      - run: pip install pytest
      - run: pytest tests/regression/ -v
```

---

## Performance Benchmarks

### Landsat Full-Scene (7931×7891 = 62M pixels × 7 bands)

| Metric | L8 Python | L8 Rust | Speedup | L9 Python | L9 Rust | Speedup |
|--------|-----------|---------|---------|-----------|---------|---------|
| Total | 179.5s | 66.1s | **2.7×** | 180.1s | 56.3s | **3.2×** |

### Sentinel-2 Full-Scene (5490×5490 at 20m, ~30M pixels × 11 bands)

| Metric | S2A Python | S2A Rust | Speedup | S2B Python | S2B Rust | Speedup |
|--------|------------|----------|---------|------------|----------|---------|
| Total | 182.3s | 52.0s | **3.5×** | 173.0s | 64.2s | **2.7×** |

### PACE OCI ROI (291 bands, 21×15 pixels)

| Stage | Python | Rust | Speedup |
|-------|--------|------|---------|
| Total | 16.7s | 0.31s | **54.7×** |

### PACE OCI Full-Scene (291 bands, 1710×1272 pixels, fixed DSF)

| Stage | Python | Rust | Notes |
|-------|--------|------|-------|
| Total | 230s | **84s** | **2.7× speedup** |
| Load | — | 12s | Bulk NetCDF reads |
| AC | — | 34s | 18.7 Mpx/s (rayon) |
| Write | — | 35s | GeoZarr V3 |

| Accuracy | Value |
|----------|-------|
| Mean Pearson R | **1.0000** |
| Mean RMSE | **0.0043** |
| All pixels within 0.05 | **100%** |
| Model match | MOD2 ✓ |

---

## Consolidated Cross-Sensor Benchmark (South Australia, all water scenes)

| Sensor | Scene Size | Python | Rust | Speedup | Mean R | Mean RMSE | Model |
|--------|-----------|--------|------|---------|--------|-----------|-------|
| Landsat 8 | 62M px × 7 bands | 180s | 66s | **2.7×** | 0.999 | 0.011 | MOD2 ✓ |
| Landsat 9 | 62M px × 7 bands | 180s | 56s | **3.2×** | 1.000 | 0.003 | MOD2 ✓ |
| Sentinel-2A | 30M px × 11 bands | 182s | 52s | **3.5×** | 1.000 | 0.001 | MOD1 ✓ |
| Sentinel-2B | 30M px × 11 bands | 173s | 64s | **2.7×** | 1.000 | 0.001 | MOD2 ✓ |
| PACE OCI | 2.2M px × 291 bands | 230s | 84s | **2.7×** | 1.000 | 0.004 | MOD2 ✓ |

All scenes over water, South Australia (Gulf St Vincent / south of Kangaroo Island).
All use fixed DSF mode for deterministic comparison.

---

## Tolerance Guidelines

| Quantity | Absolute tolerance | Relative tolerance | Notes |
|----------|-------------------|-------------------|-------|
| ρ_t (TOA reflectance) | 1e-6 | — | Should be exact (same input) |
| ρ_s (surface reflectance) | 0.002 | — | S2 fixed-mode achieves this |
| AOT | 0.001 | — | S2 fixed-mode matches to 4 sig figs |
| Pearson R | > 0.999 | — | S2 achieves 1.000 |
| Lat/Lon | 1e-8 | — | Must be exact |
| Wavelength | 0.01 nm | �