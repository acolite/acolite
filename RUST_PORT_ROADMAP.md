# ACOLITE Rust Port Roadmap

**Issue**: [#118 - Progress performance optimization work by porting to Rust](https://github.com/acolite/acolite/issues/118)

**Goal**: Port ACOLITE atmospheric correction software from Python to Rust to achieve significant performance improvements through parallelization, reduced I/O overhead, and efficient memory management.

## Background

The current Python implementation has been optimized with multiprocessing (acolite-mp branch by @woodcockr), achieving:
- **Sentinel-2**: 197s vs 791s (4x speedup)
- **Landsat**: 99s vs 312s (3.15x speedup)

However, limitations remain:
- Per-band parallelism only
- GIL contention
- Intermediate NetCDF writes
- Memory constraints with high-resolution data

## Objectives

1. **Performance**: 5-10x improvement over current Python implementation
2. **Memory Efficiency**: Keep data in RAM, use Zarr for intermediate persistence only when necessary
3. **Parallelization**: Block/chunked processing with true parallelism
4. **Compatibility**: Maintain output format compatibility with existing ACOLITE workflows
5. **Maintainability**: Clean, modular Rust architecture

---

## Phase 1: Foundation & Architecture (Weeks 1-4)

### 1.1 Project Setup
- [ ] Initialize Rust workspace with cargo
- [ ] Set up CI/CD pipeline (GitHub Actions)
- [ ] Configure testing framework
- [ ] Establish code style guidelines (rustfmt, clippy)
- [ ] Create project documentation structure

### 1.2 Core Dependencies Selection
- [ ] **Geospatial**: `gdal`, `proj`, `geo-types`
- [ ] **Array Processing**: `ndarray`, `ndarray-parallel`
- [ ] **I/O Formats**:
  - NetCDF: `netcdf` or `hdf5-rust`
  - Zarr: `zarr` or custom implementation
  - GeoTIFF: `geotiff` via GDAL bindings
- [ ] **Parallelism**: `rayon`, `tokio` (for async I/O)
- [ ] **Scientific Computing**: `nalgebra`, `approx`
- [ ] **Image Processing**: `image`, custom implementations

### 1.3 Architecture Design
- [ ] Define module structure mirroring Python codebase
- [ ] Design trait-based sensor abstraction
- [ ] Create processing pipeline architecture
- [ ] Design memory management strategy (chunking, streaming)
- [ ] Define error handling patterns
- [ ] Specify configuration system (TOML/YAML)

**Deliverables**:
- Project skeleton with build system
- Architecture documentation
- Dependency evaluation report

---

## Phase 2: Core Data Structures & I/O (Weeks 5-8)

### 2.1 Fundamental Types
- [ ] Implement `Projection` type (wrapping PROJ)
- [ ] Create `GeoTransform` structure
- [ ] Build `BandData<T>` generic array wrapper
- [ ] Implement `Metadata` dictionary structure
- [ ] Create `SensorConfig` trait and implementations

### 2.2 I/O Layer
- [ ] **NetCDF Reader/Writer**
  - Read/write with CF conventions
  - Attribute handling
  - Chunked reading for large files
- [ ] **Zarr Support**
  - Read/write Zarr arrays
  - Cloud storage backends (S3, Azure)
  - Compression (Blosc, Zstd)
- [ ] **GeoTIFF Support**
  - Read via GDAL
  - Write COG (Cloud Optimized GeoTIFF)
- [ ] **Generic Band Reader**
  - Abstraction over different formats
  - Lazy loading
  - Tile-based access

### 2.3 Sensor Metadata Parsers
Priority sensors (Landsat 8/9, Sentinel-2):
- [ ] Landsat metadata parser (MTL files)
- [ ] Sentinel-2 metadata parser (XML)
- [ ] Generic sensor trait implementation

**Deliverables**:
- Core data structures library
- I/O module with tests
- Metadata parsers for priority sensors

---

## Phase 3: Atmospheric Correction Core (Weeks 9-14)

### 3.1 Radiometric Calibration
- [ ] DN to radiance conversion
- [ ] Radiance to reflectance conversion
- [ ] Solar irradiance integration (F0 tables)
- [ ] Sun-Earth distance correction

### 3.2 Gas Correction
- [ ] Ozone absorption (k_o3 tables)
- [ ] Water vapor correction
- [ ] Gas transmittance LUT interpolation
- [ ] Ancillary data integration (OBPG)

### 3.3 Rayleigh Correction
- [ ] Rayleigh scattering LUT
- [ ] Pressure correction (DEM integration)
- [ ] Molecular atmosphere modeling

### 3.4 Aerosol Correction (DSF Algorithm)
- [ ] Dark Spectrum Fitting implementation
- [ ] Aerosol LUT interpolation
- [ ] AOT optimization
- [ ] Path reflectance calculation
- [ ] Adjacency effect correction (RAdCor)

### 3.5 LUT Management
- [ ] LUT download from acolite_luts repository
- [ ] LUT caching system
- [ ] Efficient LUT interpolation (multi-dimensional)
- [ ] Memory-mapped LUT access

**Deliverables**:
- Atmospheric correction library
- LUT management system
- Unit tests with reference data
- Performance benchmarks

---

## Phase 4: Sensor Support (Weeks 15-20)

### 4.1 Priority Sensors (Phase 4a)
- [ ] **Landsat 8/9 OLI**
  - L1 conversion
  - Band-specific processing
  - Thermal bands (TIRS)
- [ ] **Sentinel-2 MSI**
  - L1C conversion
  - Multi-resolution handling (10m/20m/60m)
  - Detector stitching

### 4.2 Secondary Sensors (Phase 4b)
- [ ] Landsat 5/7 (TM/ETM+)
- [ ] Sentinel-3 OLCI
- [ ] PlanetScope (Dove/SuperDove)
- [ ] VIIRS

### 4.3 Sensor Abstraction
- [ ] Generic `Sensor` trait
- [ ] RSR (Relative Spectral Response) handling
- [ ] Band mapping system
- [ ] Sensor-specific corrections

**Deliverables**:
- Sensor modules for priority platforms
- Sensor abstraction framework
- Validation against Python implementation

---

## Phase 5: Parallelization & Performance (Weeks 21-24)

### 5.1 Chunked Processing
- [ ] Implement tile-based processing
- [ ] Overlap handling for spatial filters
- [ ] Memory-efficient chunking strategy
- [ ] Progress tracking

### 5.2 Parallel Execution
- [ ] Band-level parallelism (rayon)
- [ ] Tile-level parallelism
- [ ] Pipeline parallelism (producer-consumer)
- [ ] Thread pool configuration

### 5.3 Memory Optimization
- [ ] Zero-copy operations where possible
- [ ] Streaming I/O for large files
- [ ] Zarr intermediate storage
- [ ] Memory profiling and optimization

### 5.4 Performance Benchmarking
- [ ] Benchmark suite against Python implementation
- [ ] Profiling (flamegraph, perf)
- [ ] Optimization of hot paths
- [ ] SIMD optimization where applicable

**Deliverables**:
- Parallel processing framework
- Performance benchmarks
- Optimization report

---

## Phase 6: Derived Products & Output (Weeks 25-28)

### 6.1 Water Quality Parameters
- [ ] Turbidity (Nechad, Dogliotti algorithms)
- [ ] Chlorophyll-a (OC, RE algorithms)
- [ ] SPM (Suspended Particulate Matter)
- [ ] CDOM absorption
- [ ] QAA (Quasi-Analytical Algorithm)

### 6.2 Output Generation
- [ ] NetCDF L2R/L2W output
- [ ] GeoTIFF export
- [ ] COG (Cloud Optimized GeoTIFF)
- [ ] RGB image generation
- [ ] Metadata preservation

### 6.3 Visualization
- [ ] Color scaling
- [ ] RGB composites
- [ ] PNG/JPEG export
- [ ] Color tables

**Deliverables**:
- Derived products library
- Output generation module
- Visualization tools

---

## Phase 7: Integration & Testing (Weeks 29-32)

### 7.1 End-to-End Pipeline
- [ ] Command-line interface
- [ ] Configuration file support
- [ ] Batch processing
- [ ] Error handling and logging

### 7.2 Validation
- [ ] Pixel-level comparison with Python ACOLITE
- [ ] Statistical validation
- [ ] Visual inspection of outputs
- [ ] Edge case testing

### 7.3 Documentation
- [ ] API documentation (rustdoc)
- [ ] User guide
- [ ] Processing examples
- [ ] Migration guide from Python

### 7.4 Python Bindings (Optional)
- [ ] PyO3 bindings for Python interop
- [ ] Python package distribution
- [ ] Hybrid workflows

**Deliverables**:
- Complete Rust ACOLITE implementation
- Comprehensive test suite
- User documentation
- Validation report

---

## Phase 8: Advanced Features (Weeks 33-36)

### 8.1 Cloud Processing
- [ ] S3/Azure Blob storage support
- [ ] Zarr cloud backends
- [ ] Distributed processing (optional)

### 8.2 Additional Sensors
- [ ] Hyperspectral sensors (PRISMA, DESIS, EMIT)
- [ ] Thermal sensors (TACT integration)
- [ ] Additional multispectral platforms

### 8.3 Advanced Corrections
- [ ] Sun glint correction
- [ ] Whitecaps correction
- [ ] Bidirectional effects
- [ ] DEM shadow masking

**Deliverables**:
- Cloud processing capabilities
- Extended sensor support
- Advanced correction algorithms

---

## Phase 9: Deployment & Maintenance (Ongoing)

### 9.1 Release Preparation
- [ ] Binary releases (Linux, macOS, Windows)
- [ ] Docker containers
- [ ] Package managers (cargo, conda)
- [ ] Version tagging and changelog

### 9.2 Community
- [ ] GitHub repository setup
- [ ] Issue templates
- [ ] Contributing guidelines
- [ ] Code of conduct

### 9.3 Maintenance
- [ ] Bug fixes
- [ ] Performance improvements
- [ ] Dependency updates
- [ ] New sensor additions

**Deliverables**:
- Production-ready releases
- Community infrastructure
- Maintenance plan

---

## Success Metrics

### Performance Targets
- **5-10x speedup** over Python implementation
- **<50% memory usage** compared to Python
- **Linear scaling** with core count up to 16 cores
- **Sub-second** processing for 10m Sentinel-2 tiles

### Quality Targets
- **<1e-6 difference** in float32 outputs vs Python
- **100% test coverage** for core algorithms
- **Zero data loss** in format conversions
- **CF-compliant** NetCDF outputs

### Usability Targets
- **Drop-in replacement** for Python ACOLITE CLI
- **<5 minute** setup time
- **Clear error messages** and logging
- **Comprehensive documentation**

---

## Risk Mitigation

### Technical Risks
1. **GDAL Rust bindings limitations**
   - Mitigation: Evaluate alternatives, contribute to gdal-rust
2. **LUT interpolation performance**
   - Mitigation: Profile early, consider GPU acceleration
3. **Memory management complexity**
   - Mitigation: Extensive testing, memory profiling tools

### Project Risks
1. **Scope creep**
   - Mitigation: Phased approach, MVP first
2. **Compatibility issues**
   - Mitigation: Continuous validation against Python
3. **Resource constraints**
   - Mitigation: Focus on priority sensors first

---

## Team & Resources

### Required Skills
- Rust systems programming
- Geospatial data processing
- Atmospheric correction algorithms
- Performance optimization
- Scientific computing

### Infrastructure
- CI/CD (GitHub Actions)
- Test data repository
- Benchmark infrastructure
- Documentation hosting

---

## Timeline Summary

| Phase | Duration | Key Deliverables |
|-------|----------|------------------|
| 1. Foundation | 4 weeks | Architecture, dependencies |
| 2. Core I/O | 4 weeks | Data structures, I/O layer |
| 3. AC Core | 6 weeks | Atmospheric correction |
| 4. Sensors | 6 weeks | Landsat, Sentinel-2 support |
| 5. Performance | 4 weeks | Parallelization, optimization |
| 6. Products | 4 weeks | Derived products, output |
| 7. Integration | 4 weeks | Testing, validation, docs |
| 8. Advanced | 4 weeks | Cloud, additional sensors |
| 9. Deployment | Ongoing | Releases, maintenance |

**Total**: ~36 weeks (9 months) for full implementation

---

## Next Steps

1. **Immediate** (Week 1):
   - Set up Rust project structure
   - Evaluate and select core dependencies
   - Create initial architecture document

2. **Short-term** (Weeks 2-4):
   - Implement basic I/O for NetCDF
   - Create Landsat metadata parser
   - Prototype DSF algorithm core

3. **Medium-term** (Weeks 5-12):
   - Complete atmospheric correction pipeline
   - Implement Landsat 8/9 support
   - Validate against Python outputs

4. **Long-term** (Weeks 13+):
   - Add Sentinel-2 support
   - Optimize performance
   - Expand sensor coverage

---

## References

- [ACOLITE Python Repository](https://github.com/acolite/acolite)
- [ACOLITE-MP Branch](https://github.com/woodcockr/acolite/blob/acolite-mp/README-ACOLITE-MP.md)
- [Issue #118](https://github.com/acolite/acolite/issues/118)
- DSF Algorithm Papers:
  - Vanhellemont & Ruddick 2018
  - Vanhellemont 2019a, 2019b
  - Vanhellemont 2020c

---

## Appendix A: Module Structure

```
acolite-rs/
├── Cargo.toml
├── src/
│   ├── lib.rs
│   ├── main.rs                 # CLI entry point
│   ├── core/
│   │   ├── mod.rs
│   │   ├── projection.rs       # Geospatial projections
│   │   ├── metadata.rs         # Metadata handling
│   │   └── band.rs             # Band data structures
│   ├── io/
│   │   ├── mod.rs
│   │   ├── netcdf.rs           # NetCDF I/O
│   │   ├── zarr.rs             # Zarr I/O
│   │   ├── geotiff.rs          # GeoTIFF I/O
│   │   └── reader.rs           # Generic reader
│   ├── sensors/
│   │   ├── mod.rs
│   │   ├── trait.rs            # Sensor trait
│   │   ├── landsat.rs          # Landsat support
│   │   ├── sentinel2.rs        # Sentinel-2 support
│   │   └── rsr.rs              # RSR handling
│   ├── ac/                     # Atmospheric correction
│   │   ├── mod.rs
│   │   ├── rayleigh.rs
│   │   ├── gas.rs
│   │   ├── aerosol.rs
│   │   ├── dsf.rs              # Dark Spectrum Fitting
│   │   └── lut.rs              # LUT management
│   ├── parameters/             # Derived products
│   │   ├── mod.rs
│   │   ├── turbidity.rs
│   │   ├── chlorophyll.rs
│   │   └── qaa.rs
│   ├── output/
│   │   ├── mod.rs
│   │   ├── netcdf.rs
│   │   ├── geotiff.rs
│   │   └── rgb.rs
│   ├── parallel/
│   │   ├── mod.rs
│   │   ├── chunking.rs
│   │   └── pipeline.rs
│   └── utils/
│       ├── mod.rs
│       ├── config.rs
│       └── logging.rs
├── tests/
│   ├── integration/
│   └── validation/
├── benches/
│   └── performance.rs
└── docs/
    ├── architecture.md
    └── user_guide.md
```

---

## Appendix B: Key Algorithms to Port

### Priority 1 (Core AC)
1. Dark Spectrum Fitting (DSF)
2. Rayleigh correction
3. Gas transmittance
4. Aerosol LUT interpolation
5. TOA to surface reflectance

### Priority 2 (Derived Products)
1. Nechad turbidity
2. Dogliotti turbidity
3. OC chlorophyll
4. Red-edge chlorophyll
5. QAA

### Priority 3 (Advanced)
1. RAdCor adjacency correction
2. Sun glint correction
3. TACT thermal processing
4. DEM shadow masking

---

**Document Version**: 1.0  
**Last Updated**: 2026-03-05  
**Author**: @whatnick  
**Status**: Draft
