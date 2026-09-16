# ACOLITE-RS Roadmap — Updated 2026-08-11

[![Rust CI](https://github.com/whatnick/acolite/actions/workflows/rust.yml/badge.svg?branch=feature/rust-port)](https://github.com/whatnick/acolite/actions/workflows/rust.yml)
[![Benchmarks](https://github.com/whatnick/acolite/actions/workflows/benchmark-suite.yml/badge.svg?branch=feature/rust-port)](https://github.com/whatnick/acolite/actions/workflows/benchmark-suite.yml)
[![S3 OLCI](https://img.shields.io/badge/S3_OLCI-44×_faster-brightgreen)](BENCHMARK_RESULTS.md)
[![S2 MSI](https://img.shields.io/badge/S2_MSI-pending-yellow)](BENCHMARK_RESULTS.md)
[![Landsat 8/9](https://img.shields.io/badge/Landsat_8/9-9×_faster-brightgreen)](BENCHMARK_RESULTS.md)

## Executive Summary

### Current State (as of 2026-08-11)
- **Rust LOC**: 13,014 across 64 `.rs` files
- **Python LOC** (target): 55,424 across 528 `.py` files
- **Ported**: 4 sensor loaders (Landsat 8/9, Sentinel-2, Sentinel-3 OLCI, PACE OCI) + full LUT-DSF atmospheric correction pipeline
- **Coverage**: ~9,758 Python LOC equivalently ported (~18% of codebase by functionality)
- **Tests**: 135 Rust + 195 Python regression = **330 total tests**
- **Elapsed time**: ~12 active days (2026-03-05 to 2026-03-17), agent-assisted with Kiro + Copilot dual-agent harness
- **Last upstream sync**: 2026-08-11 (commit `64a02ff3`)

### New Python Work Since Last Port Update (2026-07-30 → 2026-08-11)

13 commits upstream adding **82 net new lines** including:

| Change | Files | Category | Port Impact |
|--------|-------|----------|-------------|
| hDSF glint correction guard for hyperspectral | `hdsf/hdsf.py` | AC core | ✅ Ported (glint_correct_guarded) |
| nc_to_geotiff nodata handling (gdal.Warp) | `output/nc_to_geotiff.py` | Writer | ✅ Ported (COG nodata fix) |
| WorldView XML namespace-agnostic parsing | `worldview/metadata_parse.py` | Sensor loader | Phase D (WorldView) |
| WorldView PGC bundle handling improvements | `worldview/l1_convert.py` | Sensor loader | Phase D (WorldView) |
| WorldView PAN band identifiers (WV03/WV04) | `worldview/metadata_parse.py` | Sensor loader | Phase D (WorldView) |
| WorldView non-zero nodata for PGC imagery | `worldview/l1_convert.py` | Sensor loader | Phase D (WorldView) |
| WorldView multi-band file bundle support | `worldview/l1_convert.py` | Sensor loader | Phase D (WorldView) |
| WorldView half-pixel projection alignment | `worldview/l1_convert.py` | Sensor loader | Phase D (WorldView) |
| WorldView global_dims from projection data | `worldview/l1_convert.py` | Sensor loader | Phase D (WorldView) |
| Pléiades skip missing tile | `pleiades/l1_convert.py` | Sensor loader | Phase D (Pléiades) |
| Pléiades moved up components initialisation | `pleiades/bundle_test.py` | Sensor loader | Phase D (Pléiades) |
| Hyperspectral sensor check in hDSF | `hdsf/hdsf.py` | AC core | ✅ Ported |

**Impact on Rust port**:
- ✅ 2 changes ported in this sync (glint guard + COG nodata)
- 10 changes affect WorldView/Pléiades loaders (not yet ported — tracked for Phase D)
- **Phase D note**: WorldView loader must implement XML namespace-agnostic parsing (use `getElementsByTagNameNS('*', tag)` equivalent in Rust XML parser), PGC stretch detection beyond just 'mr', and support for single multi-band TIFF files as bundles

---

## Revised Total Scope Assessment

### Python LOC by Port Category

| Category | Python LOC | Rust Status | Complexity |
|----------|-----------|-------------|------------|
| **Already ported (4 sensors + DSF AC)** | 9,758 | ✅ 13,014 Rust LOC | Done |
| **Remaining sensor loaders (36 sensors)** | 13,508 | Not started | Medium (repetitive) |
| **Core orchestration** (`acolite/acolite/`) | 7,890 | Partial (pipeline.rs) | High |
| **Shared utilities** (`shared/`) | 4,780 | Partial (~30%) | Medium |
| **Adjacency correction** (RAdCor + GLAD) | 2,915 | Not started | Very High (physics) |
| **GEE integration** | 2,694 | Not started | Medium (may skip) |
| **TACT** (thermal correction) | 2,053 | Not started | High (libRadtran) |
| **API/download** | 1,808 | Partial (CMR/STAC) | Medium |
| **AerLUT** | 1,568 | ✅ Done | — |
| **AC core** | 1,987 | ✅ Done | — |
| **Output/writer** | 1,333 | Partial (COG/GeoZarr) | Medium |
| **RTM** (Hydrolight + libRadtran) | 1,037 | Not started | High |
| **DEM** | 988 | Not started | Medium |
| **GEM** (generic extract/merge) | 954 | Not started | Medium |
| **Parameters** (L2W products) | 916 | Not started | Medium |
| **Masking** | 592 | Not started | Low |
| **Glint** | 384 | ✅ Done (glint.rs) | — |
| **hDSF** | 337 | Not started | High (new) |
| **Other** (zarr, map, convert, custom) | ~700 | Partial | Low |

**Total remaining**: ~41,400 Python LOC to port (including new upstream additions)

---

## Time Estimation — Agent-Assisted Port

### Methodology
Based on observed velocity from Phase A–C (the 12-day sprint):
- **Ported**: ~9,758 Python LOC → 13,014 Rust LOC
- **Produced**: 330 regression tests
- **Effective rate**: ~815 Python LOC ported per active day (including tests + validation)
- **Rust expansion factor**: ~1.33× (Rust is more verbose + explicit error handling + tests)

### Complexity Multipliers

| Complexity | Multiplier | Reasoning |
|-----------|-----------|-----------|
| Repetitive sensor loaders | 0.6× | Template-based, agent excels at these |
| Standard utilities | 0.8× | Straightforward NumPy → ndarray |
| Core orchestration | 1.5× | State management, settings, multi-path dispatch |
| Physics algorithms (RAdCor, TACT, hDSF) | 2.0× | Numerical precision critical, external deps |
| External system integration (GEE, libRadtran) | 2.5× | May need FFI or redesign |

### Phase-by-Phase Estimate

| Phase | Scope (Python LOC) | Complexity | Agent-Days | Calendar Weeks |
|-------|-------------------|-----------|-----------|----------------|
| **D: Additional Sensor Loaders** | 13,508 | 0.6× | 10 | 2–3 |
| **E: Production Hardening** | ~3,000 | 1.0× | 4 | 1 |
| **F: L2W Water Products** | 916 + 337 (hDSF) | 1.5× | 3 | 1 |
| **G: Core Orchestration** | 7,890 | 1.5× | 15 | 3–4 |
| **H: Shared Utilities** | ~3,300 (remaining) | 0.8× | 4 | 1 |
| **I: Adjacency Correction** | 2,915 | 2.0× | 7 | 2 |
| **J: TACT + RTM** | 3,090 | 2.0× | 8 | 2 |
| **K: DEM + Masking + GEM** | 2,534 | 1.0× | 3 | 1 |
| **L: API/Download + Output** | 2,141 (remaining) | 0.8× | 3 | 1 |
| **M: GEE** (optional) | 2,694 | 2.5× | 8 | 2 |
| **N: Integration + CI + Docs** | — | 1.0× | 5 | 1–2 |
| **Total** | **~41,400** | | **70 agent-days** | **17–22 weeks** |

### Without GEE (recommended — GEE is a cloud service, better wrapped than ported)

| | Agent-Days | Calendar Weeks |
|---|-----------|----------------|
| **Total (excluding GEE)** | **62 agent-days** | **15–20 weeks** |

### Summary Estimate

| Metric | Value |
|--------|-------|
| **Already done** | 12 active agent-days (Phase A–C) |
| **Remaining (excl. GEE)** | ~62 agent-days |
| **Total project** | **~74 agent-days** |
| **Calendar time** (1 dev, 4 days/week effective) | **18–20 weeks** |
| **Calendar time** (2 agents in parallel) | **10–12 weeks** |
| **Expected final Rust LOC** | ~55,000–60,000 |
| **Expected final test count** | ~800–1,000 |

---

## Updated Phase Plan

### Phase D — Remaining Sensor Loaders (36 sensors)
**Estimate: 10 agent-days**

Batch by I/O pattern (agent can template these):

| Batch | Sensors | Pattern | Days |
|-------|---------|---------|------|
| D1: GeoTIFF-based | WorldView, Pléiades, Planet, IKONOS, Deimos, Formosat, GF, Huanjing, OpenCosmos, SDGSAT, AMAZONIA | GDAL GeoTIFF + XML metadata | 3 |
| D2: HDF5-based | PRISMA, DESIS, EnMAP, EMIT, HICO, CHRIS, HYPERION, HyperField, HYPSO, Tanager | HDF5 reader, wavelength tables | 3 |
| D3: NetCDF-based | VIIRS, GOCI, SeaHawk, EarthCare, SeaDAS | NetCDF + orbit geometry | 2 |
| D4: Geostationary | GOES, Himawari, SEVIRI, FCI | Full-disk subsetting, segment stitching | 2 |

### Phase E — Production Hardening
**Estimate: 4 agent-days**

- [ ] NetCDF L2R output (matching Python format) — **blocks Phase F**
- [ ] Ancillary spatial+temporal interpolation (MERRA2 bilinear)
- [ ] DEM-derived pressure (barometric formula from SRTM/Copernicus DEM)
- [ ] Streaming processing for large scenes (tiled I/O)
- [ ] CLI matching Python ACOLITE settings file format
- [ ] Settings TOML → ProcessingConfig mapping

### Phase F — L2W Water Products + hDSF
**Estimate: 3 agent-days**

- [ ] Port `acolite_l2w.py` parameter dispatch
- [ ] Chlorophyll-a (OC algorithms: OC3, OC4, OC5, OC6)
- [ ] Chlorophyll-a colour-ratio (`chl_crat`)
- [ ] TSS (Nechad 2010, Dogliotti 2015)
- [ ] Turbidity (Dogliotti, Nechad)
- [ ] QAA (quasi-analytical algorithm)
- [ ] MALH
- [ ] hDSF (hyperspectral DSF variant — new 336 LOC)
- [ ] Regression: derived products Rust vs Python

### Phase G — Core Orchestration
**Estimate: 15 agent-days** (most complex phase)

- [ ] `acolite_run.py` → multi-scene pipeline orchestrator
- [ ] `acolite_l1r.py` → L1 reader dispatch (identify sensor → call loader)
- [ ] `acolite_l2r.py` → L2R processing (2,244 LOC — largest single file)
- [ ] `identify_bundle.py` → input detection heuristics
- [ ] `acolite_map.py` → mapping/visualization output
- [ ] `acolite_flags.py` → quality flags
- [ ] `acolite_pans.py` → pansharpening
- [ ] `parameter_scaling.py` / `parameter_discretisation.py`
- [ ] Settings system (TOML-based, sensor defaults, user overrides)
- [ ] Logging framework

### Phase H — Shared Utilities (remaining)
**Estimate: 4 agent-days**

- [ ] Projection handling (setup, merge, sub, warp)
- [ ] NetCDF read/write/extract utilities
- [ ] Polygon/ROI/limit processing
- [ ] RSR convolution (rsr_dict, rsr_convolute)
- [ ] Sun position / geometry calculations
- [ ] WOPP (water optical properties)
- [ ] Authentication (EarthData .netrc)
- [ ] SFTP upload
- [ ] Array utilities (remaining: convolve, derivative)

### Phase I — Adjacency Correction
**Estimate: 7 agent-days** (physics-heavy)

- [ ] RAdCor algorithm (1,979 LOC — physics-based adjacency)
- [ ] GLAD algorithm (482 LOC)
- [ ] PSF modelling
- [ ] Validation against Python outputs

### Phase J — TACT + RTM
**Estimate: 8 agent-days** (external dependency complexity)

- [ ] libRadtran integration (FFI or process spawning)
- [ ] TACT thermal correction pipeline
- [ ] ERA5/ECMWF profile download and parsing
- [ ] Emissivity estimation (NDVI-based)
- [ ] Hydrolight RTM interface
- [ ] TACT simulation LUT generation

### Phase K — DEM + Masking + GEM
**Estimate: 3 agent-days**

- [ ] Copernicus DEM download + tile stitching
- [ ] SRTM DEM support
- [ ] DEM shadow masking
- [ ] Land/water mask (OCM)
- [ ] Hillshade
- [ ] GEM (scene download, extract, merge)

### Phase L — API/Download + Output (remaining)
**Estimate: 3 agent-days**

- [ ] EarthExplorer API
- [ ] CDSE (Copernicus Data Space) API
- [ ] GOCI download
- [ ] GeoTIFF output (nc_to_geotiff)
- [ ] NetCDF compression
- [ ] Reprojection/cropping utilities

### Phase M — GEE Integration (optional)
**Estimate: 8 agent-days** (recommended: wrap Python or skip)

- [ ] Google Earth Engine interface
- [ ] AGH (ACOLITE GEE Hybrid) — 1,229 LOC
- [ ] Scene finding via GEE

### Phase N — Integration, CI, Documentation
**Estimate: 5 agent-days**

- [ ] Full end-to-end pipeline integration tests
- [ ] CI/CD pipeline (GitHub Actions)
- [ ] Performance regression CI (fail on >10% throughput drop)
- [ ] Binary releases (cross-compilation)
- [ ] User documentation + migration guide
- [ ] PyO3 Python bindings (optional — drop-in replacement)

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|-----------|
| Upstream Python keeps diverging | Medium | Merge regularly (already done 2026-07-30); track diff weekly |
| Numerical precision in RAdCor/TACT | High | Early validation against Python reference scenes |
| libRadtran FFI complexity | Medium | Use process spawning (like Python does) not native FFI |
| GDAL dependency for all sensors | Low | Already feature-gated; tiff-crate fallback for simple formats |
| GEE integration infeasible in Rust | Low | Keep as Python wrapper or skip entirely |
| hDSF is still experimental upstream | Low | Port after API stabilises; only 336 LOC |

---

## Architecture (unchanged)

```
src/
├── auth/           # Secure credentials (.netrc, env vars, config)
├── loader/         # INPUT: Read satellite data (4 sensors done, 36 remaining)
├── ac/             # PROCESSING: Atmospheric correction (DSF done)
├── writer/         # OUTPUT: Write results (COG + GeoZarr done)
├── core/           # Data types (band, metadata, projection)
├── sensors/        # Sensor definitions
└── (pipeline, parallel, resample, simd)

Future additions:
├── l2w/            # Water product algorithms (Phase F)
├── adjacency/      # RAdCor + GLAD (Phase I)
├── thermal/        # TACT (Phase J)
├── dem/            # DEM handling (Phase K)
└── masking/        # Land/water/cloud masking (Phase K)
```

---

## Completed ✅ (Phases A–C, 12 active days)

- 4 sensor loaders (Landsat, S2, S3, PACE)
- Full LUT-DSF atmospheric correction (fixed + tiled)
- N-dimensional LUT interpolation
- Gas transmittance (O₃, H₂O, CO₂, O₂)
- Sun glint correction (Cox-Munk + Fresnel)
- Geographic subsetting (all 4 sensors, UTM-aware)
- Ancillary download (OBPG + MERRA2)
- COG + GeoZarr output
- Sensor-specific LUT auto-download
- 135 Rust tests + 195 Python regression tests
- Performance: 2.5–7.3× speedup across all sensors

---

## Key Metrics

| Metric | Current | End Target |
|--------|---------|-----------|
| Python LOC covered | 9,758 (18%) | 55,342 (100%) |
| Rust LOC | 13,014 | ~55,000–60,000 |
| Sensors ported | 4 / 40 | 40 / 40 |
| Test count | 330 | ~800–1,000 |
| Mean speedup | 3.4× | 3–5× (target) |
| Agent-days spent | 12 | ~74 total |
