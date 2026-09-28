---
name: regression-testing
description: "Write and run regression tests comparing Rust acolite-rs against Python ACOLITE. Use when: adding regression tests for a new sensor, validating Rust output against Python baselines, debugging numerical differences between implementations, checking physics equivalence (Pearson R, RMSE, reflectance tolerances), or setting up test tiers (synthetic, integration, real data)."
argument-hint: "Sensor name or test scenario (e.g., 'landsat surface reflectance' or 'pace band count')"
---

# ACOLITE Regression Testing

## When to Use
- Adding regression tests for a newly ported sensor
- Validating Rust output matches Python ACOLITE within tolerances
- Debugging numerical differences between Rust and Python
- Setting up tiered test suites (synthetic → integration → real data)
- Running cross-implementation performance comparisons

## Test Architecture

```
tests/
├── <sensor>_e2e.rs              # Rust-only E2E (synthetic data, always runs)
├── integration_tests.rs         # Cross-sensor pipeline tests
├── regression/
│   ├── conftest.py              # Shared fixtures: build Rust, tolerances, CLI opts
│   ├── test_<sensor>_regression.py      # Tiered Python↔Rust comparison
│   └── test_<sensor>_rust_vs_python.py  # Full AC per-pixel comparison
benches/
│   └── performance.rs           # Criterion benchmarks
```

## Test Tier Design

### Tier 1 — Synthetic (always runs, no I/O)
Validate physics invariants using synthetic data:
- Rayleigh τ decreases with wavelength (λ⁻⁴ law)
- Gas transmittance ∈ (0, 1]
- Reflectance calibration: higher sun → lower reflectance
- Parallel and sequential produce identical results

### Tier 2 — Integration (needs netCDF4/GDAL)
Validate metadata parsing and structure:
- Sensor XML/MTL parsing produces expected fields
- Output format structure (GeoZarr dirs, COG tags)
- CF-1.8 metadata conventions
- Band resolution groups match sensor spec

### Tier 3 — Real Data (`--runslow`)
Validate against actual satellite scenes:
- Band count Python = Rust (exact)
- Per-band mean reflectance within tolerance
- Spectral ordering preserved (Blue > SWIR for Rayleigh)
- Spatial extent within ROI bounds

### Tier 4 — Cross-Implementation (`--runslow`)
Full atmospheric correction comparison:
- Per-pixel Pearson R > 0.80 (typically > 0.99)
- RMSE < 0.02 for surface reflectance
- > 95% of pixels within 0.05 absolute difference
- Both select same aerosol model (Continental/Maritime)

## Procedure: Add Regression Tests for a New Sensor

### Step 1: Create Rust E2E test
Create `tests/<sensor>_e2e.rs`:
```rust
#[test]
fn test_<sensor>_band_count() { /* exact match */ }
#[test]
fn test_<sensor>_parallel_determinism() { /* two runs identical */ }
#[test]
fn test_<sensor>_reflectance_range() { /* ρ ∈ (-0.1, 2.0) */ }
```

### Step 2: Create Python regression test
Create `tests/regression/test_<sensor>_regression.py`:
```python
# Tier 1: Synthetic tests (no pytest.mark.slow)
def test_<sensor>_band_wavelengths(): ...
def test_<sensor>_calibration_linearity(): ...

# Tier 3: Real data tests
@pytest.mark.slow
def test_<sensor>_band_count_real(sensor_file): ...
@pytest.mark.slow
def test_<sensor>_reflectance_stats(sensor_file): ...
```

### Step 3: Create Rust vs Python comparison
Create `tests/regression/test_<sensor>_rust_vs_python.py`:
```python
@pytest.mark.slow
def test_<sensor>_pearson_r(): ...     # R > 0.80
@pytest.mark.slow
def test_<sensor>_rmse(): ...          # RMSE < 0.02
@pytest.mark.slow
def test_<sensor>_spectral_order(): ...
```

### Step 4: Add conftest fixtures
In `tests/regression/conftest.py`, add:
- CLI option: `--<sensor>-file` for real data path
- Environment variable: `ACOLITE_<SENSOR>_TEST_DIR`
- Fixture: `<sensor>_file` with `pytest.mark.slow` skip

### Step 5: Update benchmarks
In `benches/performance.rs`, add throughput benchmarks for the new sensor.

## Tolerance Reference

| Quantity | Tolerance | Notes |
|----------|-----------|-------|
| ρ_t (TOA reflectance) | 1e-6 | Must be exact (same input data) |
| ρ_s (surface reflectance) | 1e-4 | AC algorithm differences acceptable |
| AOT (aerosol optical thickness) | 0.02 | Optimization convergence varies |
| Band count | exact | No tolerance |
| Wavelength | 0.01 nm | Must match file metadata |
| Lat/Lon | 1e-8 | Must be exact |
| Pearson R (spatial pattern) | > 0.80 | Typically > 0.99 for validated sensors |
| RMSE (surface reflectance) | < 0.02 | Per-band pixel comparison |

## Running Tests

```bash
# Rust E2E (no external data)
cargo test --test <sensor>_e2e -- --nocapture

# Python Tier 1+2 (no real data)
pytest tests/regression/test_<sensor>_regression.py -v

# Python Tier 3+4 (real data)
pytest tests/regression/ -v --runslow --<sensor>-file /path/to/data

# Benchmarks
cargo bench

# Full regression suite
cargo test && pytest tests/regression/ -v
```

## Existing Sensor Coverage

| Sensor | Rust E2E | Python Regression | Rust vs Python | Benchmarks |
|--------|----------|-------------------|----------------|------------|
| Landsat 8/9 | ✅ 15 tests | ✅ 3 tiers | ✅ 13 tests | ✅ |
| PACE OCI | ✅ 5 tests | ✅ 4 tiers | ✅ 14 tests | — |
| Sentinel-2 | ✅ 18 tests | ✅ 3 tiers | — | ✅ |

## Reference Files
- Full strategy: [REGRESSION_TESTING_ROADMAP.md](../../REGRESSION_TESTING_ROADMAP.md)
- Port status: [RUST_PORT_ROADMAP.md](../../RUST_PORT_ROADMAP.md)
- Conftest fixtures: `tests/regression/conftest.py`
- Performance benchmarks: `benches/performance.rs`
