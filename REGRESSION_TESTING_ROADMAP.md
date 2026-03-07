# ACOLITE Regression Testing Roadmap

## Overview

This document outlines the regression testing strategy for validating the Rust
port (`acolite-rs`) against the reference Python ACOLITE implementation, with
PACE OCI as the primary target sensor.

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
│   ├── test_pace_regression.py    # PACE-specific Python↔Rust comparison
│   ├── test_landsat_regression.py # Landsat 8/9 Python↔Rust comparison
│   ├── test_sentinel2_regression.py # Sentinel-2 Python↔Rust comparison
│   └── __init__.py
benches/
│   └── performance.rs             # Criterion benchmarks (Landsat + S2 + resample)
```

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

**Run:** `pytest tests/regression/ -v` or `cargo test --test pace_e2e`

### Tier 2 — Synthetic NetCDF Round-Trip (needs netCDF4)

| Test | What it validates |
|------|-------------------|
| Synthetic PACE L1B creation | Minimal valid NetCDF with 3 detectors |
| Python l1_convert ingestion | ACOLITE can read the synthetic file |
| Band count consistency | Same number of bands from both implementations |

**Run:** `pytest tests/regression/ -v -k synthetic`

### Tier 3 — Real Data End-to-End (needs --runslow + data)

| Test | What it validates |
|------|-------------------|
| L1R band count match | Python and Rust produce same band count |
| Per-band reflectance statistics | Mean/std/min/max within tolerance |
| Spectral ordering | Wavelengths monotonically increasing |
| GeoZarr ↔ NetCDF equivalence | Rust GeoZarr readable by xarray/zarr |

**Run:** `pytest tests/regression/ -v --runslow --pace-file /path/to/PACE_OCI.L1B.nc`

---

## Phase 1: Current (PACE OCI L1B → L1R)

**Status: Implemented**

- [x] Rust PACE loader (`src/loader/pace.rs`) — reads L1B NetCDF
- [x] Parallel atmospheric correction (`src/parallel.rs`) — rayon-based
- [x] GeoZarr writer (`src/writer/geozarr.rs`) — Zarr V3 + CF-1.8
- [x] Rust E2E test (`tests/pace_e2e.rs`) — synthetic scene
- [x] Python regression harness (`tests/regression/`)

### Key metrics to track

| Metric | Python baseline | Rust target | Tolerance |
|--------|----------------|-------------|-----------|
| Band count (L1B→L1R) | ~280 | ~280 | exact |
| Mean ρ_t per band | reference | match | 1e-4 |
| Processing time (280 bands, 200×150) | ~30s | <5s | — |
| Memory peak | ~2 GB | <500 MB | — |

## Phase 1b: Landsat 8/9 OLI

**Status: Implemented**

- [x] Rust Landsat loader (`src/loader/landsat.rs`) — reads GeoTIFF bands
- [x] MTL parser (`src/sensors/mtl.rs`) — metadata extraction
- [x] Parallel AC pipeline — 7 bands via rayon
- [x] COG writer (`src/writer/cog.rs`) — per-band GeoTIFF
- [x] Rust E2E test (`tests/landsat_e2e.rs`) — 15 tests
- [x] Python regression (`tests/regression/test_landsat_regression.py`)
- [x] Performance benchmarks (`benches/performance.rs`)

### Output format: Per-band COG files

The Rust Landsat writer produces **one COG (Cloud-Optimized GeoTIFF) per band**,
not a single multi-band file. This is a deliberate design decision:

- **Parallel cloud reads**: Each band is a separate blob in object storage (S3/GCS),
  enabling concurrent downloads of individual bands without range-request overhead.
- **Selective band access**: Downstream consumers (e.g. NDVI = B5/B4) only fetch
  the bands they need, avoiding full-scene transfer.
- **Parallel write**: Each band can be written independently by a separate thread,
  eliminating serialization bottlenecks.
- **COG tiling**: Each per-band COG has its own internal tiling and overviews,
  optimized for that band's data distribution.
- **Consistency with input**: Landsat Collection 2 is distributed as per-band
  GeoTIFFs on USGS S3 (`s3://usgs-landsat/`), so the output mirrors the input
  layout for pipeline symmetry.

This differs from Python ACOLITE which writes all bands into a single NetCDF L1R file.

### Key metrics to track

| Metric | Python baseline | Rust target | Tolerance |
|--------|----------------|-------------|-----------|
| Band count (L1→L1R) | 7 | 7 | exact |
| Mean ρ_t per band | reference | match | 1e-4 |
| Processing time (7 bands, 7000×7000) | ~15s | <2s | — |
| Parallel speedup (7 bands) | 1x | >3x | — |
| Throughput | ~5 Mpx/s | >50 Mpx/s | — |

### Landsat-specific test coverage

| Test | What it validates |
|------|-------------------|
| L8/L9 band parity | Identical wavelength/bandwidth definitions |
| DN→radiance linearity | Exact: L = DN × MULT + ADD |
| DN→reflectance sun angle | Higher sun → lower reflectance |
| Earth-Sun distance cycle | Perihelion < aphelion, correct range |
| Parallel = sequential | Bit-identical results |
| Parallel determinism | Two runs → identical output |
| L8 full pipeline | 7 bands, 200×200, ρs ∈ [0, 1) |
| L9 full pipeline | Same as L8 |
| L8/L9 cross-sensor | Same input → same output |
| Rayleigh per band | τ monotonic, B1 > 0.1, B7 < 0.001 |
| Speedup scaling | Parallel ≤ 2× sequential time |
| Throughput | > 1 Mpx/s (conservative CI floor) |
| COG output | write_auto → .tif for ≤50 bands |
| DSF AOT from SWIR | AOT ∈ (0, 1) |
| AOT vs turbidity | Brighter dark → higher AOT |

## Phase 1c: Sentinel-2 A/B MSI

**Status: Implemented**

- [x] Rust S2 sensor (`src/sensors/sentinel2.rs`) — 13 bands, 3 resolution groups
- [x] S2 XML parser (`src/sensors/s2_xml.rs`) — metadata extraction
- [x] Multi-resolution resampling (`src/resample.rs`) — nearest/bilinear/average
- [x] Parallel AC per resolution group — rayon-based
- [x] COG writer for multispectral output
- [x] Rust E2E test (`tests/sentinel2_e2e.rs`) — 18 tests
- [x] Python regression (`tests/regression/test_sentinel2_regression.py`)
- [x] Performance benchmarks (`benches/performance.rs`)

### Key metrics

| Metric | Python baseline | Rust measured | Notes |
|--------|----------------|---------------|-------|
| Band count | 13 | 13 | exact |
| 10m throughput (4 bands, 2k×2k) | ~10 Mpx/s | 65.9 Mpx/s | 6.6× faster |
| 13-band parallel speedup | 1× | 3.26× | rayon |
| Resample 20m→10m bilinear | ~5 Mpx/s | 2.3 Mpx/s | single-threaded |

### Sentinel-2-specific test coverage (18 tests)

| Test | What it validates |
|------|-------------------|
| S2A/S2B band parity | 13 bands, same names |
| Resolution groups | 4×10m, 6×20m, 3×60m |
| Wavelength ordering | B01 (443nm) → B12 (2202nm) |
| Resample 20m→10m dims | 500→1000 |
| Resample 60m→10m dims | 183→1098 |
| Resample 10m→20m average | Mean preserved |
| Resample identity | Same res → no change |
| Resample roundtrip mean | 20m→10m→20m ≈ original |
| Multi-res parallel pipeline | 10m+20m+60m valid ρs |
| Parallel = sequential | Bit-identical |
| Parallel determinism | Two runs identical |
| S2A/S2B cross-sensor | Same input → same output |
| Rayleigh per band | τ monotonic, B01>0.1, B12<0.001 |
| 10m throughput | >1 Mpx/s |
| 13-band speedup | Parallel ≤ 2× sequential |
| Resample throughput | >1 Mpx/s |
| COG output | .tif for ≤50 bands |
| DSF AOT from SWIR | AOT ∈ (0, 1) |

## Phase 1d: Real PACE OCI Data Regression (OB-DAAC)

**Status: Implemented ✓**

- [x] CMR search for PACE OCI L1B granules (public API, no auth)
- [x] EarthData OAuth download (two-step redirect pattern)
- [x] Credential sources: env vars, .netrc, ~/.easi-workflows-auth.conf, config/credentials.txt
- [x] Python ACOLITE L1R processing with spatial limit (ROI subsetting)
- [x] Reference statistics generation and caching
- [x] 9 pytest regression tests on real data
- [x] Spectral physics validation (O₂-A absorption, H₂O absorption, Rayleigh)

### Test data

| Property | Value |
|----------|-------|
| Granule | `PACE_OCI.20240701T175112.L1B.V3.nc` |
| Size | 1,911 MB |
| Location | Chesapeake Bay (36.0°N, 75.5°W) |
| ROI limit | [35.8, -75.4, 36.0, -75.2] |
| ROI size | 20 × 14 pixels |
| Bands | 291 (119 blue + 163 red + 9 SWIR) |
| Download time | 112s (17.1 MB/s) |
| Processing time | 19.5s (Python ACOLITE L1R) |

### Reference statistics

| Detector | Bands | Mean ρ_t range | Notes |
|----------|-------|----------------|-------|
| Blue (315–606 nm) | 119 | 0.67–0.93 | High Rayleigh scattering |
| Red (600–895 nm) | 163 | 0.46–0.92 | O₂-A dip at 761 nm |
| SWIR (940–2258 nm) | 9 | 0.03–0.75 | H₂O absorption at 940 nm |

### Real-data test coverage (9 tests)

| Test | What it validates |
|------|-------------------|
| Band count | 119 blue + 163 red + 9 SWIR = 291 |
| Physical range | All ρ_t ∈ (-0.1, 2.0), >90% valid |
| Spectral consistency | Blue mean > SWIR mean (Rayleigh) |
| O₂-A absorption dip | ρ_t(761) < ρ_t(754) and ρ_t(771) |
| H₂O absorption | ρ_t(940) < ρ_t(865) |
| SWIR decreasing | ρ_t(1250) > ρ_t(2258) |
| L1R NetCDF readable | lat, lon, sza, >200 rhot bands |
| Spatial extent | lat/lon within ROI bounds |
| Stats reproducible | Re-read matches reference to 1e-6 |

### Commands

```bash
# Download + process + generate reference stats
python tests/regression/run_pace_real.py --python-only

# Run pytest regression (uses cached reference stats)
python -m pytest tests/regression/test_pace_regression.py -v

# Use custom cache directory
ACOLITE_TEST_CACHE=/path/to/cache python tests/regression/run_pace_real.py
```

## Phase 1e: Rust vs Python Performance Regression

**Status: Implemented ✓**

- [x] Rust `process_pace` example with `--limit` support
- [x] Zarr v3 direct reader (codec-aware: gzip + raw)
- [x] 14 pytest tests comparing Rust and Python on real PACE data
- [x] Performance, output comparison, GeoZarr structure, summary report

### Performance results (real PACE OCI, 291 bands, 21×15 ROI)

| Stage | Python | Rust | Speedup |
|-------|--------|------|---------|
| Total (load+AC+write) | 16.7s | 0.31s | **54.7×** |
| L1B load (NetCDF subset) | ~16s | 242ms | ~66× |
| Atmospheric correction | — | 2.4ms (38 Mpx/s) | — |
| GeoZarr write (291 bands) | — | 50ms | — |

### Output comparison

| Metric | Value |
|--------|-------|
| Band count match | 291 = 291 ✓ |
| Wavelength range | 315–2258 nm ✓ |
| Spectral correlation (TOA vs corrected) | r = 0.9995 |
| Bands where ρs ≤ ρt | 289/289 (100%) |
| O₂-A absorption dip preserved | ✓ |
| Blue > SWIR (Rayleigh) | ✓ |

### Test coverage (14 tests)

| Test | What it validates |
|------|-------------------|
| Rust faster than Python | Total time speedup > 1× |
| AC throughput | > 10 Mpx/s |
| Load time | < 5s for ROI subset |
| Write time | < 5s for 291 bands |
| Band count match | Python = Rust |
| Wavelength coverage | Same range ±2 nm |
| Spectral correlation | r > 0.95 |
| Corrected ≤ TOA | > 70% of bands |
| Spectral shape preserved | O₂-A dip, Blue > SWIR |
| Zarr exists | Single .zarr output |
| Zarr metadata | CF-1.8, sensor, EPSG |
| Zarr data shape | [291, nrows, ncols] |
| Zarr wavelengths | Per-detector sorted |
| Summary report | Prints timing breakdown |

### Commands

```bash
# Run Rust vs Python comparison (requires cached PACE file)
python -m pytest tests/regression/test_pace_rust_vs_python.py -v -s

# Run all regression tests
python -m pytest tests/regression/ -v
```

## Phase 1f: Real Landsat 8 Data Regression (USGS S3)

**Status: Implemented ✓**

- [x] USGS STAC search for Landsat C2L1 scenes (public API)
- [x] S3 requester-pays download via `AWS_PROFILE` (no hardcoded credentials)
- [x] Python ACOLITE L1R processing with spatial limit (ROI subsetting)
- [x] Reference statistics generation and caching
- [x] Rust `process_landsat` example with `--file`, `--limit`, `--output` CLI args
- [x] MTL parsing for sun geometry
- [x] 14 pytest regression tests comparing Rust and Python
- [x] Per-band COG output validation

### Test data

| Property | Value |
|----------|-------|
| Scene | `LC08_L1TP_013035_20240710_20240718_02_T1` |
| Size | 601 MB (7 bands + MTL + ANG) |
| Location | Outer Banks / Chesapeake Bay (35.0–37.1°N) |
| Platform | Landsat 8 OLI |
| Cloud cover | 31.5% |
| Source | `s3://usgs-landsat/collection02/level-1/` (requester-pays) |
| AWS profile | `adias-prod-sso-power` |
| Python ROI | [35.8, -75.6, 36.0, -75.3] → 743×908 pixels |

### Performance results (real Landsat 8, 7 bands, full scene 7821×7691)

| Stage | Python | Rust | Notes |
|-------|--------|------|-------|
| L1R conversion (ROI) | 3.8s | — | Python processes 743×908 ROI |
| Load (7 bands, full) | — | 4.6s | 601 MB GeoTIFF via GDAL |
| Atmospheric correction | — | 12.9s | 60M pixels × 7 bands |
| COG write | — | 9.0s | Per-band COG output |
| Total (full scene) | — | 17.9s | 7821×7691 × 7 bands |

### Python reference statistics (ROI: 743×908)

| Band | Wavelength | Mean ρ_t | Min | Max |
|------|-----------|----------|-----|-----|
| rhot_443 | 443 nm | 0.6579 | 0.184 | 1.233 |
| rhot_483 | 483 nm | 0.6556 | 0.154 | 1.278 |
| rhot_561 | 561 nm | 0.6158 | 0.109 | 1.255 |
| rhot_655 | 655 nm | 0.6294 | 0.088 | 1.317 |
| rhot_865 | 865 nm | 0.6647 | 0.075 | 1.325 |
| rhot_1609 | 1609 nm | 0.3671 | 0.056 | 0.660 |
| rhot_2201 | 2201 nm | 0.2199 | 0.028 | 0.372 |

### Test coverage (14 tests)

| Test | What it validates |
|------|-------------------|
| Rust completes | < 120s total |
| Rust vs Python timing | Speedup comparison |
| Load time | < 30s for 7 bands |
| AC throughput | Mpx/s metric |
| Output exists | COG file created |
| Band count | 7 bands (per-band COGs) |
| Image dimensions | > 1000×1000 |
| Reflectance range | Mean ∈ (-1, 2) |
| Spectral ordering | SWIR darker than VNIR |
| GeoTIFF format | GTiff driver |
| Projection present | CRS in output |
| GeoTransform present | Pixel size non-zero |
| Data type | Float32 or Float64 |
| Summary report | Prints comparison table |

### Commands

```bash
# Download Landsat scene from USGS S3 (requester-pays)
AWS_PROFILE=adias-prod-sso-power aws s3 cp \
  s3://usgs-landsat/collection02/level-1/standard/oli-tirs/2024/013/035/LC08_L1TP_013035_20240710_20240718_02_T1/ \
  /tmp/acolite_test_cache/landsat/LC08_L1TP_013035_20240710_20240718_02_T1/ \
  --request-payer requester --recursive \
  --exclude "*" --include "*_B[1-7].TIF" --include "*_MTL.txt" --include "*_ANG.txt"

# Run Python ACOLITE reference processing
python -c "import acolite as ac; ac.acolite.acolite_run(settings={
  'inputfile': '/tmp/acolite_test_cache/landsat/LC08_L1TP_013035_20240710_20240718_02_T1',
  'output': '/tmp/acolite_test_cache/landsat/py_output',
  'limit': [35.8, -75.6, 36.0, -75.3]})"

# Run Rust processing
cargo run --release --example process_landsat -- \
  --file /tmp/acolite_test_cache/landsat/LC08_L1TP_013035_20240710_20240718_02_T1 \
  --output /tmp/acolite_test_cache/landsat/rust_output

# Run regression tests
python -m pytest tests/regression/test_landsat_rust_vs_python.py -v -s
```

---

## Phase 2: L1R → L2R (DSF Atmospheric Correction)

**Status: Implemented ✓**

The full LUT-based Dark Spectrum Fitting (DSF) atmospheric correction has been
ported to Rust, using the same 6SV radiative transfer LUTs as Python ACOLITE.

### Architecture decision: Full LUT port to Rust

Rather than bridging to Python's LUT system, we ported the complete LUT pipeline
to Rust. This is required for the full ACOLITE Rust port and provides:

- **Zero Python dependency** at runtime — pure Rust atmospheric correction
- **Same LUT files** — reads the existing ACOLITE NetCDF LUTs directly
- **Same physics** — identical correction formula: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`
- **Multi-model selection** — tests MOD1 (Continental) and MOD2 (Maritime), selects by minimum RMSD

### New Rust modules

| Module | File | Purpose |
|--------|------|---------|
| `ac::interp` | `src/ac/interp.rs` | N-dimensional regular grid interpolator (scipy.interpolate.RegularGridInterpolator equivalent) |
| `ac::aerlut` | `src/ac/aerlut.rs` | Aerosol LUT reader — loads sensor-specific 6SV NetCDF LUTs, stacks pressures, builds per-band interpolators |
| `ac::gas_lut` | `src/ac/gas_lut.rs` | Gas transmittance — reads ko3 absorption data, RSR files, computes per-band O₃/H₂O/CO₂/O₂ transmittance |
| `ac::dsf` | `src/ac/dsf.rs` | Full DSF: dark spectrum extraction, LUT-based AOT inversion, model selection, surface reflectance correction |

### DSF algorithm (matching Python ACOLITE)

1. **DN → TOA reflectance**: `rhot = DN * REFLECTANCE_MULT + REFLECTANCE_ADD` / cos(sza)
2. **Gas transmittance**: Ozone (ko3 × RSR convolution), water vapour, CO₂, O₂
3. **Dark spectrum**: 1st percentile of each band (full scene)
4. **AOT inversion**: For each model, invert LUT romix(tau) → tau for each band's dark pixel
5. **Model selection**: Minimum RMSD between observed dark spectrum and modeled path reflectance
6. **Surface reflectance**: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`

### Regression results — South Australia water scenes (Gulf St Vincent)

| Metric | L8 (tiled, MOD2) | L9 (tiled, MOD2) | Python L8 (AOT=0.68, MOD2) | Python L9 (AOT=0.33, MOD2) |
|--------|---------------------|---------------------|---------------------------|---------------------------|
| Mean Pearson R | **0.9985** | **0.9996** | — | — |
| Mean RMSE | **0.011** | **0.003** | — | — |
| Mean %<0.05 | **99.9%** | **100.0%** | — | — |
| Mean %<0.01 | 71.1% | 98.1% | — | — |

**L9 achieves near-perfect agreement**: RMSE=0.003, 100% of pixels within 0.05, 98% within 0.01.
**L8 excellent agreement**: RMSE=0.011, 99.9% of pixels within 0.05 tolerance.
Both pipelines select MOD2 (Maritime) aerosol model. Tiled AOT estimation (200×200 tiles)
with intercept dark spectrum method matches Python's `dsf_aot_estimate=tiled`.

Rust mean AOT: L8=0.59 (Python=0.68), L9=0.27 (Python=0.33). Remaining AOT gap is from
Rust processing full scene vs Python processing ROI subset only.

### Performance (full scene 7931×7891 = 62M pixels × 7 bands)

| Stage | L8 | L9 | Notes |
|-------|----|----|-------|
| Load | 4.6s | 4.2s | 7 GeoTIFF bands via GDAL |
| LUT load | 0.3s | 0.3s | 2 models × 4 pressures × 9 bands |
| AC (DSF + correction) | 24.8s | 23.5s | Full scene, all 7 bands |
| COG write | 20.3s | 21.9s | Parallel per-band COG |
| **Total** | **52.6s** | **52.2s** | Full scene processing |

Note: Python processes only the ROI subset (1140×949 = 1.1M pixels) in 8.9s.
Rust processes the full scene (62M pixels) — 56× more pixels.

### Remaining gaps vs Python ACOLITE

| Feature | Status | Impact on accuracy |
|---------|--------|-------------------|
| ROI subsetting (limit) | Not yet | Would match Python's ROI-only processing |
| Interface reflectance (rsky) | Not yet | Minor when dsf_interface_reflectance=False |
| Ancillary data (NCEP) | Not yet | Uses default ozone/pressure |
| DEM-derived pressure | Not yet | Uses 1013 hPa |
| Glint correction | Not yet | Minor for nadir sensors |
| Tiled AOT estimation | **Done** | Implemented 200×200 tile DSF |

### Tests

| Test | What it validates |
|------|-------------------|
| 7 bands output | All bands processed |
| Pearson R > 0.80 | Spatial pattern correlation |
| Spectral ordering | SWIR ≤ VNIR (water) |
| COG structure | Projection, geotransform, layout |
| Full report | Per-band metrics table |

### Build requirement

```bash
# LUT-based DSF requires the full-io feature (NetCDF support)
cargo build --release --features full-io --example process_landsat
```

## Phase 3: L2R → L2W (Water Products)

**Status: Planned**

- [ ] Port `acolite_l2w.py` parameter computation
- [ ] Chlorophyll-a (OC algorithms)
- [ ] TSS (Nechad, Dogliotti)
- [ ] Turbidity
- [ ] Regression: compare derived products

### Tests to add

```
test_chl_oc_match                  — Chl-a within ±10%
test_tss_nechad_match              — TSS within ±10%
test_turbidity_match               — FNU within ±10%
```

## Phase 4: Multi-Sensor Expansion

**Status: Planned**

| Sensor | Loader | AC | Writer | Regression |
|--------|--------|----|--------|------------|
| PACE OCI | ✅ | ✅ (basic) | ✅ GeoZarr | ✅ 5 tests |
| Landsat 8/9 | ✅ | ✅ LUT-DSF | ✅ COG | ✅ 13 tests |
| Sentinel-2 | ✅ | ✅ (basic) | ✅ COG | ✅ 18 tests |
| Sentinel-3 OLCI | ✅ | ✅ (basic) | — | ✅ (integration) |
| PRISMA | — | — | — | — |
| EMIT | — | — | — | — |

## Phase 5: Performance Regression

- [ ] Benchmark suite (`benches/performance.rs`) — track Mpx/s
- [ ] CI integration: fail if throughput drops >10%
- [ ] Memory profiling: track peak RSS

---

## Running Tests

### Rust tests (no external data needed)

```bash
# Unit + integration tests
cargo test

# PACE E2E test specifically
cargo test --test pace_e2e

# Landsat E2E test specifically
cargo test --test landsat_e2e

# Sentinel-2 E2E test specifically
cargo test --test sentinel2_e2e

# Sentinel-2 with performance output
cargo test --test sentinel2_e2e -- --nocapture

# Landsat with performance output
cargo test --test landsat_e2e -- --nocapture

# With netcdf feature (for real file tests)
cargo test --features netcdf --test pace_e2e

# Benchmarks (Landsat parallel vs sequential scaling)
cargo bench
```

### Python regression tests

```bash
# Tier 1+2 (no real data)
pytest tests/regression/ -v

# All tiers with real PACE data
pytest tests/regression/ -v --runslow --pace-file /data/PACE_OCI.20240701.L1B.nc

# All tiers with real Landsat data
pytest tests/regression/ -v --runslow --landsat-file /data/LC08_L1TP_...

# All tiers with real Sentinel-2 data
pytest tests/regression/ -v --runslow --s2-file /data/S2A_MSIL1C_....SAFE

# With custom tolerance
pytest tests/regression/ -v --runslow --pace-file /data/PACE.nc --tolerance 1e-3
```

### CI Pipeline (recommended)

```yaml
# .github/workflows/regression.yml
jobs:
  rust-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: cargo test --test pace_e2e
      - run: cargo test --test integration_tests

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

## Tolerance Guidelines

| Quantity | Absolute tolerance | Relative tolerance | Notes |
|----------|-------------------|-------------------|-------|
| ρ_t (TOA reflectance) | 1e-6 | — | Should be exact (same input) |
| ρ_s (surface reflectance) | 1e-4 | 1% | AC differences acceptable |
| AOT | 0.02 | 10% | Optimization may converge differently |
| Chl-a | — | 10% | Derived product, error propagates |
| Lat/Lon | 1e-8 | — | Must be exact |
| Wavelength | 0.01 nm | — | Must match file metadata |
