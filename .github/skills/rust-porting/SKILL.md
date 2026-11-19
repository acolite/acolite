---
name: rust-porting
description: "Port Python ACOLITE modules to Rust acolite-rs. Use when: porting a sensor loader, implementing atmospheric correction, translating NumPy to ndarray, creating Rust equivalents of Python functions, adding new sensor support, or implementing LUT-based algorithms. Covers DN→reflectance calibration, gas transmittance, Rayleigh correction, DSF inversion, and output writing."
argument-hint: "Name the Python module or function to port (e.g., 'acolite/landsat/l1_convert.py' or 'gas_transmittance')"
---

# Python-to-Rust Porting for ACOLITE

## When to Use
- Porting a Python sensor loader to Rust (e.g., `acolite/prisma/` → `src/loader/prisma.rs`)
- Translating atmospheric correction algorithms from Python to Rust
- Adding a new sensor definition in `src/sensors/`
- Implementing LUT interpolation or radiative transfer functions in Rust
- Converting NumPy array operations to `ndarray` equivalents

## Architecture Mapping

```
Python module              →  Rust location
─────────────────────────────────────────────
acolite/<sensor>/          →  src/loader/<sensor>.rs + src/sensors/<sensor>.rs
acolite/ac/rayleigh.py     →  src/ac/rayleigh.rs
acolite/ac/gas_transmittance.py → src/ac/gas.rs + src/ac/gas_lut.rs
acolite/ac/optimise_aot_*  →  src/ac/dsf.rs
acolite/aerlut/            →  src/ac/aerlut.rs
acolite/acolite/acolite_l1r.py → src/pipeline.rs (L1R stage)
acolite/acolite/acolite_l2r.py → src/pipeline.rs (L2R stage)
acolite/output/            →  src/writer/
```

## Porting Procedure

### Step 1: Understand the Python code
1. Read the target Python module completely
2. Identify all NumPy operations and their Rust equivalents
3. Note any global state, settings dict access, or file I/O
4. List all dependencies (LUT files, config, ancillary data)

### Step 2: Design the Rust interface
1. Define error types using `thiserror` (never `unwrap()` in library code)
2. Use `ndarray::Array2<f32>` where Python uses 2D float arrays
3. Feature-gate NetCDF/HDF5 I/O: `#[cfg(feature = "full-io")]`
4. Prefer `rayon` for parallelism — ensure deterministic output

### Step 3: Implement with physics equivalence
Key formulas that MUST match Python exactly:
- **TOA reflectance**: `rhot = (DN * MULT + ADD) / cos(sza)`
- **Gas correction**: `rhot_gc = rhot / tt_gas` (per-band transmittance)
- **Surface reflectance**: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`
- **Rayleigh**: Use the same 6SV LUT files — do NOT reimplement the scattering model

### Step 4: Write tests
Add tests in the appropriate location:
- Unit tests: `#[cfg(test)] mod tests` inside the module
- E2E tests: `tests/<sensor>_e2e.rs`
- Regression: `tests/regression/test_<sensor>_regression.py`

## Common Python → Rust Translations

| Python (NumPy) | Rust (ndarray) |
|----------------|----------------|
| `np.zeros((h, w))` | `Array2::<f32>::zeros((h, w))` |
| `arr[arr < 0] = 0` | `arr.mapv_inplace(\|v\| v.max(0.0))` |
| `np.percentile(arr, 1)` | Sort + index (no built-in) |
| `np.interp(x, xp, fp)` | `src/ac/interp.rs` regular grid interpolator |
| `scipy.interpolate.RegularGridInterpolator` | `src/ac/interp.rs::RegularGridInterpolator` |
| `settings['key']` | `ProcessingConfig` struct fields |
| `netCDF4.Dataset(path)` | `netcdf::open(path)` behind `#[cfg(feature = "full-io")]` |

## Band Naming Conventions
- TOA: `rhot_443` (wavelength in nm)
- Surface: `rhos_443`
- Band count determines output format: >50 → GeoZarr, ≤50 → per-band COG

## Sensor Loader Template

A new sensor loader needs:
1. **Sensor definition** (`src/sensors/<sensor>.rs`): Band wavelengths, bandwidths, resolution groups
2. **Metadata parser** (`src/sensors/<sensor>_meta.rs`): Extract sun/view geometry, timestamps
3. **Loader** (`src/loader/<sensor>.rs`): Read input data → `Vec<BandData>`
4. **Register** in `src/sensors/mod.rs` and `src/loader/mod.rs`

## Reference Files
- Roadmap: [RUST_PORT_ROADMAP.md](../../RUST_PORT_ROADMAP.md) — sensor priority and status
- Regression roadmap: [REGRESSION_TESTING_ROADMAP.md](../../REGRESSION_TESTING_ROADMAP.md) — test tiers and tolerances
- Rust public API: `src/lib.rs` — all exported modules and types
- Python entry: `acolite/acolite/acolite_run.py` — main processing flow
- LUT data: `data/LUT/` — shared NetCDF LUT files read by both implementations
