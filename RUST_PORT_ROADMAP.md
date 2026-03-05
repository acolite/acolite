# ACOLITE-RS Roadmap

## Architecture

```
src/
├── auth/           # Secure credentials (.netrc, env vars, config)
│   └── credentials.rs
├── loader/         # INPUT: Read satellite data
│   ├── geotiff.rs  # GDAL-based GeoTIFF reader
│   ├── landsat.rs  # Landsat L1 scene loader
│   └── source/     # Remote data access
│       ├── cmr.rs      # NASA CMR search
│       ├── stac.rs     # STAC search
│       └── download.rs # Download (EarthData + AWS S3)
├── ac/             # PROCESSING: Atmospheric correction
│   ├── pipeline.rs → (in src/pipeline.rs)
│   ├── calibration.rs  # DN → TOA
│   ├── rayleigh.rs     # Rayleigh scattering
│   ├── gas.rs          # Gas transmittance
│   ├── dsf.rs          # Dark Spectrum Fitting
│   └── lut.rs          # LUT management
├── writer/         # OUTPUT: Write results
│   └── cog.rs      # Cloud-Optimized GeoTIFF
├── core/           # Data types
├── sensors/        # Sensor definitions
└── (pipeline, parallel, resample, simd)
```

## Current State (Post-Refactor)

- **2,026 lines** source (down from 4,233 — 52% reduction)
- **5 examples** (down from 19 — 74% reduction)
- **25 tests** passing
- **Clean architecture**: loader → ac → writer
- **Secure credentials**: zeroize-on-drop, .netrc support, no hardcoded profiles
- **No code duplication**: single search, single download, single GeoTIFF reader
- **Real data validated**: Landsat 8 scene processed (7691×7801, 18.6 Mpx/s)

## Completed ✅

### Phase A: Consolidation
- [x] Restructured into loader/ac/writer architecture
- [x] Merged 3 search implementations → `loader/source/cmr.rs`
- [x] Merged 3 download implementations → `loader/source/download.rs`
- [x] Deleted duplicate GeoTIFF reader (tiff crate version)
- [x] Deleted 4 PACE-specific IO files
- [x] Secure credential handling (zeroize, .netrc, env vars)
- [x] Consolidated 19 examples → 5
- [x] Deleted intermediate documentation

### Earlier Phases
- [x] Core data types (BandData, Metadata, Projection, GeoTransform)
- [x] Sensor definitions (Landsat 8/9, Sentinel-2, Sentinel-3, PACE)
- [x] GDAL GeoTIFF reading
- [x] COG writing
- [x] STAC + CMR search
- [x] Parallel band processing (rayon)
- [x] Real Landsat scene processed via AWS S3

## Next: Phase B — Real Atmospheric Correction

### B.1 LUT Loading (1 week)
Port LUT loading from Python ACOLITE's `acolite_luts` repo:
- [ ] Download NetCDF LUTs from https://github.com/acolite/acolite_luts
- [ ] Parse 6SV Rayleigh LUT (wavelength × SZA × VZA × RAA × pressure)
- [ ] Parse aerosol LUTs (per model: Continental, Maritime, Urban)
- [ ] Implement N-dimensional interpolation

**Python reference**: `acolite/aerlut/import_lut.py`, `acolite/aerlut/ilut.py`

### B.2 Rayleigh Correction (3 days)
Replace placeholder with real implementation:
- [ ] `ray_tau()` — Rayleigh optical thickness from LUT
- [ ] `ray_phase()` — Phase function
- [ ] `ray_refl()` — Rayleigh reflectance (path + sky)
- [ ] Pressure correction

**Python reference**: `acolite/ac/rayleigh.py` (99 lines)

### B.3 Gas Transmittance (3 days)
Replace placeholder with real implementation:
- [ ] O3 transmittance from absorption coefficients
- [ ] H2O transmittance from LUT
- [ ] O2 transmittance
- [ ] Combined gas correction

**Python reference**: `acolite/ac/gas_transmittance.py` (68 lines)

### B.4 DSF Algorithm (1 week)
Port the core atmospheric correction:
- [ ] Tile-based dark spectrum extraction
- [ ] AOT optimization per tile (minimize cost function)
- [ ] Aerosol model selection
- [ ] Path reflectance subtraction
- [ ] Diffuse transmittance correction

**Python reference**: `acolite/acolite/acolite_l2r.py` (2218 lines — core section ~500 lines)

### B.5 Calibration (2 days)
- [ ] Parse MTL for calibration coefficients
- [ ] DN → radiance → TOA reflectance chain
- [ ] Earth-Sun distance correction per scene date

**Python reference**: `acolite/landsat/read_toa.py`, `acolite/landsat/metadata_bands.py`

## Phase C — Validation (1-2 weeks)

- [ ] Process same Landsat scene with Python and Rust
- [ ] Pixel-level RMSE comparison (target: <0.001 ρ)
- [ ] Band-by-band statistical comparison
- [ ] Performance benchmarking on real data
- [ ] Edge cases: cloud-heavy, coastal, inland water

## Phase D — Python Sync Framework (Ongoing)

### Strategy
1. Pin to Python ACOLITE release tags (not main branch)
2. Track changes in critical files only:
   - `acolite/acolite/acolite_l2r.py` (core AC)
   - `acolite/ac/rayleigh.py` (Rayleigh)
   - `acolite/ac/gas_transmittance.py` (gas)
   - `acolite/aerlut/ilut.py` (LUT interpolation)
   - `acolite/landsat/l1_convert.py` (Landsat loader)
3. Automated comparison tests on tagged releases

### Python → Rust Function Mapping

| Python | Rust |
|--------|------|
| `landsat.l1_convert()` | `loader::landsat::load_landsat_scene()` |
| `acolite_l1r()` | `pipeline::Pipeline::process_band()` (calibration) |
| `acolite_l2r()` | `pipeline::Pipeline::process_band()` (AC) |
| `ac.rayleigh.ray_refl()` | `ac::rayleigh::rayleigh_correction()` |
| `ac.gas_transmittance()` | `ac::gas::gas_correction()` |
| `aerlut.ilut()` | `ac::lut::LutManager` |
| `output.nc_write()` | `writer::cog::write_cog()` |

### What NOT to Port
- GUI (`acolite_gui.py`)
- GEE integration
- TACT (thermal) — separate concern
- RAdCor (adjacency) — separate concern
- 30+ sensor `l1_convert` — only Landsat/S2/PACE initially

## Phase E — Production Hardening (2-4 weeks)

- [ ] Remove all `unwrap()` from library code
- [ ] Streaming processing (process while downloading)
- [ ] CLI matching Python ACOLITE settings files
- [ ] NetCDF L2 output matching Python format
- [ ] Sentinel-2 loader (JP2 via GDAL)
- [ ] PACE loader (HDF5/NetCDF)
- [ ] Memory optimization for large scenes
- [ ] rustdoc documentation

## Dependencies

```toml
gdal = "0.17"        # GeoTIFF read/write
ndarray = "0.15"     # Array processing
rayon = "1.8"        # Parallelism
tiff = "0.9"         # Fallback TIFF writing
reqwest = "0.11"     # HTTP
serde = "1.0"        # Serialization
chrono = "0.4"       # DateTime
thiserror = "1.0"    # Error types
zeroize = "1.8"      # Credential security
log = "0.4"          # Logging
dirs = "5.0"         # Home directory
```
