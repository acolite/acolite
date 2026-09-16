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
Python: acolite/sentinel3/  →  Rust: src/loader/sentinel3.rs, src/sensors/sentinel3.rs
```

Pipeline: `loader → ac (calibration → gas → rayleigh → DSF → glint) → writer`

## Key Files

- Roadmaps: [RUST_PORT_ROADMAP.md](../RUST_PORT_ROADMAP.md), [REGRESSION_TESTING_ROADMAP.md](../REGRESSION_TESTING_ROADMAP.md)
- Rust entry: `src/lib.rs` (public API), `src/main.rs` (CLI)
- Python entry: `acolite/acolite/acolite_run.py`
- Tests: `tests/` (Rust E2E), `tests/regression/` (Python↔Rust)
- Config: `config/defaults.txt` (Python settings)
- Glint correction: `src/ac/glint.rs` (Cox-Munk + Fresnel, `compute_rsky`, `glint_correct`)
- Multi-agent harness: `tools/agent_harness.py` (Kiro ↔ Copilot ACP orchestrator)
- Context tool: `tools/kiro_context.py` (compact resume prompt for Kiro context resets)

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

- **Physics equivalence**: Rust must produce numerically equivalent results to Python (Pearson R > 0.999, RMSE < 0.002 for S2/PACE surface reflectance, < 0.02 for Landsat)
- **Same LUT files**: Rust reads the same ACOLITE NetCDF LUTs — do not create new LUT formats
- **Correction formula**: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`
- **Band naming**: `rhot_443` (TOA), `rhos_443` (surface) — wavelength in nm
- **Tolerances**: ρ_t exact (1e-6), ρ_s within 0.002 (S2) / 0.02 (Landsat), AOT within 0.001 (S2)
- **Pearson R target**: R > 0.999 for all sensors — S3 OLCI currently at R=0.983 (AOT gap under investigation)

## Sensor Status

| Sensor | Status | Accuracy |
|--------|--------|----------|
| Landsat 8/9 | ✅ Full pipeline | R>0.998, RMSE<0.02 |
| Sentinel-2 A/B | ✅ Full pipeline | R=1.000, RMSE<0.002 |
| PACE OCI | ✅ Full pipeline | R=0.9995 |
| Sentinel-3 OLCI | ✅ Full pipeline | R=0.983 (**AOT gap: 0.076 Rust vs 0.057 Python — open issue**) |
