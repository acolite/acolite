# ACOLITE Project Guidelines

## Overview

ACOLITE is a satellite atmospheric correction toolkit. The codebase has two implementations:
- **Python ACOLITE** (`acolite/`) — Reference implementation, 48+ sensors
- **Rust acolite-rs** (`src/`) — Port in progress, targeting performance parity with physics equivalence

## Architecture

```
Python: acolite/ac/         →  Rust: src/ac/
Python: acolite/acolite/    →  Rust: src/pipeline.rs, src/parallel.rs
Python: acolite/landsat/    →  Rust: src/loader/landsat.rs, src/sensors/landsat.rs
Python: acolite/pace/       →  Rust: src/loader/pace.rs, src/sensors/pace.rs
Python: acolite/sentinel2/  →  Rust: src/sensors/sentinel2.rs, src/resample.rs
```

Pipeline: `loader → ac (calibration → gas → rayleigh → DSF) → writer`

## Key Files

- Roadmaps: [RUST_PORT_ROADMAP.md](../RUST_PORT_ROADMAP.md), [REGRESSION_TESTING_ROADMAP.md](../REGRESSION_TESTING_ROADMAP.md)
- Rust entry: `src/lib.rs` (public API), `src/main.rs` (CLI)
- Python entry: `acolite/acolite/acolite_run.py`
- Tests: `tests/` (Rust E2E), `tests/regression/` (Python↔Rust)
- Config: `config/defaults.txt` (Python settings)

## Code Style

### Rust
- Use `thiserror` for error types — no `unwrap()` in library code
- Feature-gate I/O-heavy deps: `#[cfg(feature = "full-io")]` for NetCDF
- Parallel processing via `rayon` — ensure deterministic output
- Output: per-band COG (≤50 bands) or GeoZarr (>50 bands)

### Python
- Follow existing ACOLITE conventions — no type hints unless already present
- Settings dict pattern: `settings = {'inputfile': ..., 'output': ..., 'limit': [S, W, N, E]}`

## Build & Test

```bash
# Rust
cargo test                           # All unit + E2E tests
cargo test --features full-io        # With NetCDF support
cargo bench                          # Performance benchmarks
cargo build --release --features full-io  # Release build

# Python regression
pytest tests/regression/ -v          # Tier 1+2 (no real data)
pytest tests/regression/ -v --runslow  # All tiers (needs real data)
```

## Conventions

- **Physics equivalence**: Rust must produce numerically equivalent results to Python (Pearson R > 0.99, RMSE < 0.02 for surface reflectance)
- **Same LUT files**: Rust reads the same ACOLITE NetCDF LUTs — do not create new LUT formats
- **Correction formula**: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`
- **Band naming**: `rhot_443` (TOA), `rhos_443` (surface) — wavelength in nm
- **Tolerances**: ρ_t exact (1e-6), ρ_s within 1e-4, AOT within 0.02
