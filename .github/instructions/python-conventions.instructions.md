---
description: "Use when editing Python ACOLITE files or regression tests. Covers the settings dict pattern, sensor module structure, and pytest regression test conventions."
applyTo: ["acolite/**/*.py", "tests/regression/**/*.py"]
---

# Python ACOLITE Conventions

## Settings Pattern
ACOLITE uses a flat settings dict — not classes or dataclasses:
```python
settings = {
    'inputfile': '/path/to/scene',
    'output': '/path/to/output',
    'limit': [S, W, N, E],  # ROI bounding box
    'l2w_parameters': ['chl_oc3', 'spm_nechad'],
}
ac.acolite.acolite_run(settings=settings)
```

## Regression Tests
- Mark real-data tests with `@pytest.mark.slow`
- Use conftest fixtures: `rust_binary`, `tolerance`, `<sensor>_file`
- Cache expensive outputs in `ACOLITE_TEST_CACHE` (default `/tmp/acolite_test_cache/`)
- Tiers: synthetic (always) → integration (netCDF4) → real data (`--runslow`)
- Comparison metrics: Pearson R, RMSE, MAE, % within tolerance

## Module Structure
Each sensor lives in `acolite/<sensor>/` with:
- `l1_convert.py` — Read raw data → L1R NetCDF
- `metadata.py` or inline metadata parsing
- RSR files in `data/RSR/<SENSOR>/`
