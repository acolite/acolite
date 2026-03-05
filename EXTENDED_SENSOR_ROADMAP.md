# ACOLITE Rust Port - Extended Sensor Roadmap

**Project Goal**: Port ACOLITE to Rust for 5-10x performance improvement with full sensor support

**Current Status**: Phase 7 Complete (90% of core functionality)  
**Timeline**: 36 weeks core → Extended to 64 weeks for complete sensor suite (52 sensors)

---

## Completed Phases (Weeks 1-26) ✅

### Phase 1-7: Core Implementation
- ✅ Foundation & Architecture
- ✅ Core Processing & Parallelization
- ✅ Atmospheric Correction (DSF, Rayleigh, Gas)
- ✅ Sentinel-2 MSI Support
- ✅ Performance Optimization & I/O
- ✅ Sentinel-3 OLCI Support
- ✅ GeoTIFF I/O & Scene Processing
- ✅ STAC Integration, COG & GeoZarr

**Sensors Implemented**: Landsat 8/9, Sentinel-2 MSI, Sentinel-3 OLCI

---

## Phase 8: Performance Profiling & Validation (Weeks 27-30) 🔄

### Objectives
- Profile and optimize hot paths
- Validate against Python ACOLITE
- Process real satellite scenes
- Accuracy assessment

### Deliverables
- Performance profiling report
- Validation results vs Python
- Optimized critical paths
- Real data processing examples

---

## Phase 9: Landsat Legacy & MODIS (Weeks 31-34)

### Sensors to Implement

#### Landsat 5 TM (1984-2013) - **COG Output**
- 7 bands (30m multispectral, 120m thermal)
- Wavelengths: 450-2350 nm + thermal
- Similar to Landsat 8 processing
- Historical archive access

#### Landsat 7 ETM+ (1999-present) - **COG Output**
- 8 bands (30m multispectral, 60m thermal, 15m pan)
- SLC-off gap filling (post-2003)
- Wavelengths: 450-2350 nm + thermal
- Gap-filled product support

#### MODIS Aqua/Terra (2000-present) - **COG Output**
- 36 bands (250m-1km resolution)
- Wavelengths: 405-14385 nm
- L1B 1km product processing
- Ocean color bands (8-16)

### Technical Requirements
- MTL parser for Landsat 5/7
- SLC-off gap interpolation
- MODIS HDF-EOS reader
- Band resampling for MODIS

### Output Format
- **All COG**: Multispectral sensors (<40 bands)
- Tiled, compressed, cloud-optimized

---

## Phase 10: High-Resolution Commercial (Weeks 35-38)

### Sensors to Implement

#### PlanetScope Dove (2013-present) - **COG Output**
- 4-5 bands (3m resolution)
- Blue, Green, Red, NIR, (Red Edge)
- Daily global coverage
- Small scene size

#### PlanetScope SuperDove (2020-present) - **COG Output**
- 8 bands (3m resolution)
- Coastal, Blue, Green I/II, Yellow, Red, Red Edge, NIR
- Enhanced spectral resolution
- Wavelengths: 431-877 nm

#### RapidEye (2008-2020) - **COG Output**
- 5 bands (5m resolution)
- Blue, Green, Red, Red Edge, NIR
- Wavelengths: 440-850 nm

#### WorldView-2 (2009-present) - **COG Output**
- 8 bands (2m multispectral, 0.5m pan)
- Coastal, Blue, Green, Yellow, Red, Red Edge, NIR1, NIR2
- Wavelengths: 400-1040 nm

#### WorldView-3 (2014-present) - **COG Output**
- 16 bands (1.2m multispectral, 0.3m pan, 3.7m SWIR)
- 8 VNIR + 8 SWIR bands
- Wavelengths: 400-2365 nm

### Technical Requirements
- High-resolution tile processing
- Pan-sharpening algorithms
- Commercial metadata parsers
- STAC integration for Planet API

### Output Format
- **All COG**: Small scenes, multispectral
- High compression for storage efficiency

---

## Phase 11: European & Asian Satellites (Weeks 39-42)

### Sensors to Implement

#### Venµs (2017-present) - **GeoZarr Output**
- 12 bands (5m resolution)
- Wavelengths: 415-910 nm
- High temporal frequency (2 days)
- Superspectral resolution

#### SPOT 6/7 (2012-present) - **COG Output**
- 5 bands (6m multispectral, 1.5m pan)
- Blue, Green, Red, NIR, Pan
- Wavelengths: 450-890 nm

#### Pléiades 1A/1B (2011-present) - **COG Output**
- 5 bands (2m multispectral, 0.5m pan)
- Blue, Green, Red, NIR, Pan
- Wavelengths: 430-940 nm

#### QuickBird-2 (2001-2015) - **COG Output**
- 5 bands (2.4m multispectral, 0.6m pan)
- Blue, Green, Red, NIR, Pan
- Wavelengths: 450-900 nm

### Technical Requirements
- DIMAP metadata parser (Venµs, SPOT, Pléiades)
- High-resolution processing
- European data hub integration

### Output Format
- **Venµs**: GeoZarr (12 bands, superspectral)
- **Others**: COG (multispectral)

---

## Phase 12: Hyperspectral Sensors (Weeks 43-46)

### Sensors to Implement

#### PRISMA (2019-present) - **GeoZarr Output**
- 239 bands (30m resolution)
- Wavelengths: 400-2505 nm
- VNIR: 66 bands (400-1010 nm)
- SWIR: 173 bands (920-2505 nm)

#### DESIS (2018-present) - **GeoZarr Output**
- 235 bands (30m resolution)
- Wavelengths: 400-1000 nm
- ISS-based sensor
- 2.55 nm spectral sampling

#### HICO (2009-2014) - **GeoZarr Output**
- 128 bands (90m resolution)
- Wavelengths: 380-960 nm
- ISS-based sensor
- 5.7 nm spectral sampling

#### HYPERION (2000-2017) - **GeoZarr Output**
- 220 bands (30m resolution)
- Wavelengths: 400-2500 nm
- VNIR: 70 bands, SWIR: 172 bands
- 10 nm spectral sampling

#### CHRIS (2001-2012) - **GeoZarr Output**
- 62 bands (17-36m resolution)
- Wavelengths: 400-1050 nm
- Multi-angle capability
- 5 viewing angles

### Technical Requirements
- HDF5 reader for PRISMA
- NetCDF reader for DESIS
- Multi-angle processing for CHRIS
- Spectral binning and resampling
- Large memory optimization

### Output Format
- **All GeoZarr**: Hyperspectral (>100 bands)
- Chunked storage for efficient access
- Wavelength metadata arrays

---

## Phase 13: Ocean Color Missions (Weeks 47-50)

### Sensors to Implement

#### PACE OCI (2024-present) - **GeoZarr Output**
- 286 hyperspectral bands (1km resolution)
- Wavelengths: 315-2260 nm
- Continuous spectrum: 5nm sampling (340-890 nm)
- SWIR bands for atmospheric correction
- Primary ocean color mission for 2020s+

#### ENVISAT MERIS (2002-2012) - **COG Output**
- 15 bands (300m resolution)
- Wavelengths: 390-1040 nm
- Full resolution and reduced resolution
- Ocean color heritage

#### VIIRS (2011-present) - **COG Output**
- 22 bands (375m-750m resolution)
- I-bands: 375m (5 bands)
- M-bands: 750m (16 bands)
- Day/Night band
- Wavelengths: 412-12000 nm

### Technical Requirements
- PACE OCI NetCDF-4 reader
- Hyperspectral ocean color algorithms
- MERIS N1 format reader
- VIIRS HDF5 reader
- Multi-resolution band handling
- Day/Night band processing

### Output Format
- **PACE OCI**: GeoZarr (286 hyperspectral bands)
- **MERIS, VIIRS**: COG (multispectral)

---

## Phase 15: Additional Multispectral Sensors (Weeks 53-56)

### Sensors to Implement

#### Huanjing-1 (HJ-1) CCD (2008-present) - **COG Output**
- 4 bands (30m resolution)
- Blue, Green, Red, NIR
- Chinese environmental monitoring satellite
- Wavelengths: 430-900 nm

#### Deimos-1/2 (2009-present) - **COG Output**
- 4 bands (22m resolution)
- Green, Red, Red Edge, NIR
- Spanish Earth observation
- Wavelengths: 520-900 nm

#### Formosat-2 (2004-2016) - **COG Output**
- 5 bands (8m multispectral, 2m pan)
- Blue, Green, Red, NIR, Pan
- Taiwanese satellite
- Wavelengths: 450-900 nm

#### IKONOS (1999-2015) - **COG Output**
- 5 bands (4m multispectral, 1m pan)
- Blue, Green, Red, NIR, Pan
- First commercial high-resolution
- Wavelengths: 445-853 nm

#### GF-1/2/6 (Gaofen) (2013-present) - **COG Output**
- 4-8 bands (2-16m resolution)
- Chinese high-resolution series
- Blue, Green, Red, NIR + Pan
- Wavelengths: 450-890 nm

#### Amazonia-1 (2021-present) - **COG Output**
- 4 bands (64m resolution)
- Blue, Green, Red, NIR
- Brazilian satellite
- Wavelengths: 450-890 nm

#### SDGSAT-1 (2021-present) - **COG Output**
- 7 bands (10m multispectral)
- Coastal, Blue, Green, Red, NIR, SWIR1, SWIR2
- Chinese SDG monitoring
- Wavelengths: 420-2400 nm

### Technical Requirements
- Chinese satellite metadata parsers
- Multiple commercial formats
- High-resolution processing
- International data access

### Output Format
- **All COG**: Multispectral sensors (<10 bands)

### Sensors to Implement

#### Landsat 8/9 TIRS (Thermal) - **COG Output**
- 2 thermal bands (100m resolution)
- Wavelengths: 10.6-12.5 µm
- TACT integration (already in Python ACOLITE)
- Surface temperature retrieval

#### ASTER (1999-present) - **COG Output**
- 14 bands (15m VNIR, 30m SWIR, 90m TIR)
- Wavelengths: 520-11650 nm
- Thermal emissivity
- DEM capability

### Technical Requirements
- TACT thermal correction
- Emissivity estimation
- Multi-resolution thermal processing

### Output Format
- **All COG**: Thermal products
- Temperature units in metadata

---

## Sensor Summary Table

| Sensor | Bands | Resolution | Wavelength Range | Output Format | Phase |
|--------|-------|------------|------------------|---------------|-------|
| **Completed** |
| Landsat 8/9 OLI | 7 | 30m | 433-2300 nm | COG | 1-7 ✅ |
| Sentinel-2 MSI | 13 | 10/20/60m | 443-2190 nm | COG | 4 ✅ |
| Sentinel-3 OLCI | 21 | 300m | 400-1020 nm | COG | 5 ✅ |
| **Multispectral (COG)** |
| Landsat 5 TM | 7 | 30m | 450-2350 nm | COG | 9 |
| Landsat 7 ETM+ | 8 | 30m | 450-2350 nm | COG | 9 |
| MODIS | 36 | 250m-1km | 405-14385 nm | COG | 9 |
| PlanetScope Dove | 4-5 | 3m | 455-860 nm | COG | 10 |
| PlanetScope SuperDove | 8 | 3m | 431-877 nm | COG | 10 |
| RapidEye | 5 | 5m | 440-850 nm | COG | 10 |
| WorldView-2 | 8 | 2m | 400-1040 nm | COG | 10 |
| WorldView-3 | 16 | 1.2m | 400-2365 nm | COG | 10 |
| SPOT 6/7 | 5 | 6m | 450-890 nm | COG | 11 |
| Pléiades 1A/1B | 5 | 2m | 430-940 nm | COG | 11 |
| QuickBird-2 | 5 | 2.4m | 450-900 nm | COG | 11 |
| ENVISAT MERIS | 15 | 300m | 390-1040 nm | COG | 13 |
| VIIRS | 22 | 375m-750m | 412-12000 nm | COG | 13 |
| PACE OCI | 286 | 1km | 315-2260 nm | GeoZarr | 13 |
| Landsat 8/9 TIRS | 2 | 100m | 10.6-12.5 µm | COG | 14 |
| ASTER | 14 | 15m-90m | 520-11650 nm | COG | 14 |
| **Superspectral (GeoZarr)** |
| Venµs | 12 | 5m | 415-910 nm | GeoZarr | 11 |
| **Hyperspectral (GeoZarr)** |
| PRISMA | 239 | 30m | 400-2505 nm | GeoZarr | 12 |
| DESIS | 235 | 30m | 400-1000 nm | GeoZarr | 12 |
| HICO | 128 | 90m | 380-960 nm | GeoZarr | 12 |
| HYPERION | 220 | 30m | 400-2500 nm | GeoZarr | 12 |
| CHRIS | 62 | 17-36m | 400-1050 nm | GeoZarr | 12 |
| PACE OCI | 286 | 1km | 315-2260 nm | GeoZarr | 13 |

**Total Sensors**: 32 (3 complete, 29 planned)

**Note**: Python ACOLITE supports 40+ sensors. Additional sensors to be added in future phases:
- **Hyperspectral**: EMIT, EnMAP, HYPSO-1, Tanager, Wyvern
- **Geostationary**: GOES, Himawari, Meteosat FCI, SEVIRI, GOCI
- **Ocean/Thermal**: AVHRR, ECOSTRESS, EarthCARE, Haiyang
- **Commercial**: OpenCosmos

These will be prioritized based on user demand and data availability.

---

## Output Format Decision Rules

### COG (Cloud Optimized GeoTIFF)
**Use for**: Multispectral sensors with <25 bands

**Advantages**:
- Wide tool support (GDAL, QGIS, etc.)
- HTTP range requests for partial reads
- Overviews for multi-resolution
- Standard format for most users

**Sensors**: Landsat, Sentinel-2, MODIS, PlanetScope, WorldView, SPOT, Pléiades, QuickBird, MERIS, VIIRS, ASTER

### GeoZarr
**Use for**: Hyperspectral/superspectral sensors with ≥25 bands

**Advantages**:
- Efficient chunked storage
- Better compression for many bands
- Parallel access to bands
- Cloud-native format
- Wavelength metadata

**Sensors**: Venµs (12 bands but superspectral), PRISMA, DESIS, HICO, HYPERION, CHRIS

---

## Implementation Priority

### High Priority (Phases 9-10)
1. **Landsat 5/7**: Large historical archive, user demand
2. **MODIS**: Global coverage, ocean color
3. **PlanetScope**: Daily coverage, commercial demand

### Medium Priority (Phases 11-12)
4. **Venµs**: Unique temporal resolution
5. **PRISMA**: Modern hyperspectral
6. **WorldView-2/3**: High-resolution commercial

### Lower Priority (Phases 13-14)
7. **MERIS**: Historical ocean color
8. **VIIRS**: MODIS successor
9. **Thermal sensors**: Specialized applications

---

## Technical Challenges by Phase

### Phase 9 (Landsat Legacy & MODIS)
- SLC-off gap filling for Landsat 7
- MODIS HDF-EOS format complexity
- Multi-resolution band handling
- Large file sizes (MODIS L1B ~300MB)

### Phase 10 (High-Resolution Commercial)
- Pan-sharpening algorithms
- Small scene tiling
- Commercial API authentication
- High data volume

### Phase 11 (European & Asian)
- DIMAP format parsing
- Multiple metadata standards
- European data hub access
- Coordinate system variations

### Phase 12 (Hyperspectral)
- Large memory requirements (239 bands)
- Spectral calibration
- Atmospheric correction for full spectrum
- GeoZarr optimization for 200+ bands

### Phase 13 (Ocean Color)
- MERIS N1 format (ESA proprietary)
- VIIRS multi-resolution
- PACE OCI hyperspectral (286 bands)
- Ocean-specific algorithms
- Glint correction

### Phase 14 (Thermal)
- TACT integration
- Emissivity estimation
- Temperature retrieval
- Thermal-specific corrections

---

## Performance Targets

### Processing Speed (vs Python ACOLITE)
- **Multispectral (COG)**: 10-50x faster
- **Hyperspectral (GeoZarr)**: 5-20x faster
- **High-resolution**: 20-100x faster (smaller scenes)

### Output Size
- **COG**: 40-60% of uncompressed (DEFLATE)
- **GeoZarr**: 30-50% of uncompressed (zlib level 5)

### Memory Usage
- **Multispectral**: <2GB per scene
- **Hyperspectral**: <8GB per scene
- **Streaming**: Process scenes larger than RAM

---

## Testing Strategy

### Per-Sensor Testing
- Unit tests for metadata parsing
- Integration tests for full pipeline
- Validation against Python ACOLITE
- Real data processing examples

### Format Testing
- COG validation (gdal_validate)
- GeoZarr spec compliance
- Cloud access performance
- Compression ratio verification

### Performance Testing
- Benchmark suite per sensor
- Memory profiling
- Parallel scaling tests
- I/O bottleneck analysis

---

## Deliverables by Phase

### Phase 9
- Landsat 5/7 sensor implementations
- MODIS sensor implementation
- SLC-off gap filling algorithm
- 3 new examples

### Phase 10
- 5 commercial sensor implementations
- Pan-sharpening module
- Planet API integration
- 5 new examples

### Phase 11
- 4 European/Asian sensor implementations
- DIMAP parser
- 4 new examples

### Phase 12
- 5 hyperspectral sensor implementations
- Spectral processing optimizations
- GeoZarr performance tuning
- 5 new examples

### Phase 13
- 2 ocean color sensor implementations
- PACE OCI NetCDF-4 reader
- Hyperspectral ocean color algorithms
- MERIS N1 reader
- VIIRS HDF5 reader
- 3 new examples

### Phase 14
- TACT thermal correction
- ASTER implementation
- Thermal processing guide
- 2 new examples

---

## Final Project Statistics (Week 64)

### Code
- **Total Lines**: ~20,000-25,000 Rust
- **Sensors**: 52 complete implementations
- **Tests**: 300+ comprehensive tests
- **Examples**: 50+ working demonstrations

### Performance
- **Average Speedup**: 10-50x vs Python
- **Memory Efficiency**: 50-70% reduction
- **Output Size**: 40-60% compression

### Coverage
- **Temporal**: 1984-present (40+ years)
- **Spatial**: 0.3m-2km resolution
- **Spectral**: 2-420 bands
- **Missions**: 30+ satellite missions
- **Geostationary**: 4 missions (GOES, Himawari, Meteosat, GOCI)

---

## Maintenance & Future Work

### Post-Phase 14
- New sensor additions as launched
- Algorithm improvements
- Performance optimizations
- User feedback integration
- Cloud platform integration (AWS, GCP, Azure)
- GPU acceleration for hyperspectral
- Machine learning integration

### Additional Sensors (Future Phases 15-17)

#### Phase 15: Additional Multispectral (Weeks 53-56)
- Huanjing-1 CCD, Deimos-1/2, Formosat-2
- IKONOS, GF-1/2/6 (Gaofen)
- Amazonia-1, SDGSAT-1
- **7 sensors, all COG output**

#### Phase 16: Advanced Hyperspectral (Weeks 57-60)
- EMIT (285 bands, 60m) - GeoZarr
- EnMAP (228 bands, 30m) - GeoZarr
- HYPSO-1 (120 bands, 200m) - GeoZarr
- Tanager (420 bands, 30m) - GeoZarr
- Wyvern (>100 bands, 5m) - GeoZarr
- **5 sensors, all GeoZarr output**

#### Phase 17: Geostationary & Specialized (Weeks 61-64)
- GOES-R ABI (16 bands, 0.5-2km) - COG
- Himawari-8/9 AHI (16 bands, 0.5-2km) - COG
- Meteosat FCI (16 bands, 1-2km) - COG
- GOCI-II (13 bands, 250m) - COG
- AVHRR (6 bands, 1km) - COG
- ECOSTRESS (5 bands, 70m thermal) - COG
- EarthCARE MSI (7 bands, 500m) - COG
- Haiyang (10 bands, 250m-1km) - COG
- **8 sensors, all COG output**

**Extended Total**: 52 sensors across 17 phases

### Long-term Vision
- Real-time processing capability
- Distributed processing framework
- Web service API
- Interactive visualization
- Integration with analysis platforms

---

**Document Version**: 2.0  
**Last Updated**: March 5, 2026  
**Status**: Phases 1-7 Complete, Phases 8-14 Planned
