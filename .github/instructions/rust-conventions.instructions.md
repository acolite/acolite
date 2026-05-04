---
description: "Use when editing Rust source files in the acolite-rs atmospheric correction port. Covers error handling with thiserror, feature gating for NetCDF, ndarray patterns, rayon parallelism, and physics equivalence requirements."
applyTo: "src/**/*.rs"
---

# Rust ACOLITE-RS Conventions

## Error Handling
- Use `thiserror` derive macros — never `unwrap()` or `expect()` in library code
- Propagate errors with `?` — define variants in `src/error.rs`
- `unwrap()` is acceptable only in `#[cfg(test)]` blocks and examples

## Feature Gating
- NetCDF I/O: `#[cfg(feature = "full-io")]`
- Keep core AC math available without heavy I/O deps

## Parallelism
- Use `rayon` for band-parallel processing
- Output must be deterministic (same input → identical output on re-runs)
- Test with `process_bands_parallel` vs `process_bands_sequential` parity check

## Array Operations
- Use `ndarray::Array2<f32>` for 2D raster data
- `f32` pipeline for pre-calibrated sensors (PACE), `f64` only when Python reference uses it
- SIMD helpers in `src/simd.rs` for hot loops

## Output Format
- ≤50 bands → per-band COG via `write_cog()`
- >50 bands → GeoZarr V3 via `write_geozarr()`
- `write_auto()` dispatches based on band count

## Physics
- Same LUT files as Python — read from `data/LUT/`
- Correction formula: `rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))`
- Tolerances: ρ_t exact (1e-6), ρ_s within 1e-4, AOT within 0.02
