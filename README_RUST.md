# ACOLITE-RS

Rust implementation of ACOLITE atmospheric correction for aquatic remote sensing.

## Status

🚧 **Early Development** - Phase 1 (Foundation & Architecture)

## Features

- ✅ Core data structures (Projection, Metadata, BandData)
- ✅ Sensor abstraction trait
- ✅ Landsat 8/9 sensor definitions
- ✅ Basic Rayleigh correction framework
- 🚧 NetCDF I/O (placeholder)
- 🚧 Full atmospheric correction pipeline

## Building

```bash
cargo build --release
```

## Running

```bash
cargo run --release
```

## Testing

```bash
cargo test
```

## Benchmarking

```bash
cargo bench
```

## Architecture

See [RUST_PORT_ROADMAP.md](../RUST_PORT_ROADMAP.md) for detailed implementation plan.

## License

GPL-3.0 (same as ACOLITE)
