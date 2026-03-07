# ACOLITE-RS Roadmap

## Architecture

```
src/
├── auth/           # Secure credentials (.netrc, env vars, config)
│   └── credentials.rs
├── loader/         # INPUT: Read satellite data
│   ├── geotiff.rs  # GDAL-based GeoTIFF reader
│   ├── landsat.rs  # Landsat L1 scene loader
│   ├── pace.rs     # PACE OCI L1B NetCDF reader
│   └── source/     # Remote data access
│       ├── cmr.rs      # NASA CMR search (OBDAAC + LP DAAC)
│       ├── stac.rs     # STAC search
│       └── download.rs # Download (EarthData + AWS S3)
├── ac/             # PROCESSING: Atmospheric correction
│   ├── calibration.rs, rayleigh.rs, gas.rs, dsf.rs, lut.rs
├── writer/         # OUTPUT: Write results
│   ├── mod.rs      # Auto-dispatch: >50 bands → GeoZarr, ≤50 → COG
│   ├── cog.rs      # Cloud-Optimized GeoTIFF (multispectral/superspectral)
│   └── geozarr.rs  # GeoZarr V3 (hyperspectral)
├── core/           # Data types
├── sensors/        # Sensor definitions
└── (pipeline, parallel, resample, simd)
```

## Output Format Strategy

`write_auto()` selects format based on band count:

| Category | Bands | Format | Rationale |
|----------|-------|--------|-----------|
| Hyperspectral | >50 | GeoZarr | 3D chunked array, efficient for 100s of bands |
| Superspectral | 16–50 | per-band COG | Widely compatible, manageable file count |
| Multispectral | ≤15 | per-band COG | Standard GIS workflow |

## Sensor Audit — Band Count Classification

### Hyperspectral (>50 bands) → GeoZarr

| Sensor | Bands | Status |
|--------|-------|--------|
| Tanager | 420 | Not started |
| PACE OCI | 286 | ✅ Loader + search done |
| EMIT | 285 | Not started |
| HYPERION | 242 | Not started |
| PRISMA | 239 | Not started |
| DESIS | 235 | Not started |
| EnMAP | 224 | Not started |
| HyperField | 150 | Not started |
| HICO | 128 | Not started |
| HYPSO | 120 | Not started |
| CHRIS | 62 | Not started |

### Superspectral (16–50 bands) → per-band COG

| Sensor | Bands | Status |
|--------|-------|--------|
| Aqua/Terra MODIS | 36 | Not started |
| WorldView-3 | 29 | Not started |
| VIIRS (NPP/J1/J2) | 22 | Not started |
| Sentinel-3 OLCI | 21 | Sensor def exists |
| AMAZONIA-1 WFI | 18 | Not started |
| GOES ABI | 16 | Not started |
| Himawari AHI | 16 | Not started |
| MTG-I FCI | 16 | Not started |

### Multispectral (≤15 bands) → per-band COG

| Sensor | Bands | Status |
|--------|-------|--------|
| ENVISAT MERIS | 15 | Not started |
| Sentinel-2 MSI | 13 | Sensor def exists |
| GOCI-2 | 12 | Not started |
| SEVIRI | 12 | Not started |
| Sentinel-3 SLSTR | 11 | Not started |
| Landsat 8/9 OLI | 9 | ✅ Full pipeline |
| Landsat 7 ETM+ | 8 | Not started |
| PlanetScope SD8 | 8 | Not started |
| Landsat 5 TM | 7 | Not started |
| WorldView-2 | 6–8 | Not started |
| Pléiades | 5–7 | Not started |
| QuickBird-2 | 5 | Not started |

### Priority for Rust Port

1. **Landsat 8/9** — ✅ Done
2. **PACE OCI** — ✅ Loader + search + GeoZarr writer
3. **Sentinel-2 MSI** — Sensor def exists, needs JP2 loader
4. **Sentinel-3 OLCI** — Sensor def exists, needs NetCDF loader
5. **PRISMA/DESIS/EnMAP** — Share HDF5 loader pattern
6. **EMIT** — NetCDF, similar to PACE

## Regression Testing

Full regression strategy is documented in [REGRESSION_TESTING_ROADMAP.md](REGRESSION_TESTING_ROADMAP.md).

### Rust E2E Tests

| Test file | Sensor | Tests | Command |
|-----------|--------|-------|---------|
| [tests/landsat_e2e.rs](tests/landsat_e2e.rs) | Landsat 8/9 | 15 | `cargo test --test landsat_e2e` |
| [tests/pace_e2e.rs](tests/pace_e2e.rs) | PACE OCI | 5 | `cargo test --test pace_e2e` |
| [tests/sentinel2_e2e.rs](tests/sentinel2_e2e.rs) | Sentinel-2 | 18 | `cargo test --test sentinel2_e2e` |
| [tests/integration_tests.rs](tests/integration_tests.rs) | Cross-sensor | 7 | `cargo test --test integration_tests` |
| [benches/performance.rs](benches/performance.rs) | Landsat+S2 | — | `cargo bench` |

### Python ↔ Rust Regression Tests

| Test file | Sensor | Tier | Command |
|-----------|--------|------|---------|
| [tests/regression/test_landsat_regression.py](tests/regression/test_landsat_regression.py) | Landsat 8/9 | 1-3 | `pytest tests/regression/test_landsat_regression.py -v` |
| [tests/regression/test_pace_regression.py](tests/regression/test_pace_regression.py) | PACE OCI | 1-4 | `pytest tests/regression/test_pace_regression.py -v` |
| [tests/regression/test_sentinel2_regression.py](tests/regression/test_sentinel2_regression.py) | Sentinel-2 | 1-3 | `pytest tests/regression/test_sentinel2_regression.py -v` |
| [tests/regression/test_landsat_rust_vs_python.py](tests/regression/test_landsat_rust_vs_python.py) | Landsat 8/9 | Full AC | `pytest tests/regression/test_landsat_rust_vs_python.py -v -s` |
| [tests/regression/test_pace_rust_vs_python.py](tests/regression/test_pace_rust_vs_python.py) | PACE OCI | Full AC | `pytest tests/regression/test_pace_rust_vs_python.py -v -s` |
| [tests/regression/conftest.py](tests/regression/conftest.py) | — | Fixtures | Shared pytest config, tolerances, CLI options |

### Test Tiers

- **Tier 1 — Synthetic** (always runs): Physics invariants, monotonicity, range checks
- **Tier 2 — Integration** (needs netCDF4): Metadata parsing, structure validation
- **Tier 3 — Real data** (`--runslow`): Per-band reflectance stats, band count parity
- **Tier 4 — Cross-implementation** (`--runslow`): Rust vs Python per-pixel correlation (R > 0.80, RMSE < 0.02)

### Quick Commands

```bash
# All Rust tests (no external data)
cargo test

# All Python regression (Tier 1+2, no real data)
pytest tests/regression/ -v

# Full regression with real data
pytest tests/regression/ -v --runslow \
  --pace-file /path/to/PACE_OCI.L1B.nc \
  --landsat-file /path/to/LC08_L1TP_... \
  --s2-file /path/to/S2A_MSIL1C_....SAFE
```

## Current State

- **25 tests** passing
- **5 examples** (Landsat, PACE, Sentinel-2, Sentinel-3, AWS Landsat)
- **Clean architecture**: loader → ac → writer
- **Dual output**: GeoZarr (hyperspectral) + COG (multi/superspectral)
- **Real data validated**: Landsat 8 scene + PACE OBDAAC search

## Completed ✅

### Phase A: Architecture
- [x] Loader → AC → Writer restructure
- [x] Secure credentials (zeroize, .netrc, env vars)
- [x] Single search/download/reader implementations
- [x] PACE OCI L1B NetCDF reader
- [x] OBDAAC CMR search (collection_id based)
- [x] GeoZarr writer (zarrs crate, Zarr V3, gzip, CF metadata)
- [x] Auto-dispatch writer (>50 bands → GeoZarr, ≤50 → COG)
- [x] f32 pipeline path for pre-calibrated sensors (PACE)

## Next: Phase B — Real Atmospheric Correction

### B.1 LUT Loading (1 week)
- [ ] Download NetCDF LUTs from https://github.com/acolite/acolite_luts
- [ ] Parse 6SV Rayleigh LUT (wavelength × SZA × VZA × RAA × pressure)
- [ ] Parse aerosol LUTs (Continental, Maritime, Urban)
- [ ] N-dimensional interpolation

### B.2 Rayleigh Correction (3 days)
- [ ] Replace placeholder with LUT-based implementation
- [ ] Pressure correction

### B.3 Gas Transmittance (3 days)
- [ ] O3/H2O/O2 transmittance from LUTs
- [ ] Replace all-zero k_o3 values

### B.4 DSF Algorithm (1 week)
- [ ] Tile-based dark spectrum extraction
- [ ] AOT optimization per tile
- [ ] Aerosol model selection

### B.5 Calibration (2 days)
- [ ] DN → radiance → TOA reflectance from MTL

## Phase C — Validation (1-2 weeks)
- [x] Pixel-level RMSE comparison with Python ACOLITE — see [test_landsat_rust_vs_python.py](tests/regression/test_landsat_rust_vs_python.py), [test_pace_rust_vs_python.py](tests/regression/test_pace_rust_vs_python.py)
- [x] Band-by-band statistical comparison — Landsat R>0.998, PACE r=0.9995
- [ ] Sentinel-2 Rust vs Python full AC comparison (pending LUT-DSF for S2)

## Phase D — Additional Loaders
- [ ] Sentinel-2 (JP2 via GDAL)
- [ ] Sentinel-3 OLCI (NetCDF)
- [ ] PRISMA/DESIS/EnMAP (HDF5)
- [ ] EMIT (NetCDF)

## Phase E — Production Hardening
- [ ] Remove all `unwrap()` from library code
- [ ] Streaming processing
- [ ] CLI matching Python ACOLITE settings files
- [ ] NetCDF L2 output matching Python format

## Dependencies

```toml
gdal = "0.17"        # GeoTIFF read/write
ndarray = "0.15"     # Array processing
rayon = "1.8"        # Parallelism
zarrs = "0.23"       # Zarr V3 (GeoZarr output)
tiff = "0.9"         # Fallback TIFF writing
reqwest = "0.11"     # HTTP
serde = "1.0"        # Serialization
chrono = "0.4"       # DateTime
thiserror = "1.0"    # Error types
zeroize = "1.8"      # Credential security
netcdf = "0.9"       # PACE/S3 NetCDF (optional)
```
