"""
PACE OCI Rust vs Python regression test.

Compares acolite-rs and Python ACOLITE processing of real PACE L1B data:
  - Per-band reflectance statistics
  - Processing performance (load, AC, write)
  - Output structure validation

Requires:
  - Cached PACE L1B file (run: python tests/regression/run_pace_real.py --python-only)
  - Rust binary built with --features netcdf
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CACHE_DIR = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache"))
PACE_FILE = CACHE_DIR / "PACE_OCI.20240701T175112.L1B.V3.nc"
LIMIT = [35.8, -75.4, 36.0, -75.2]
LIMIT_STR = ",".join(str(x) for x in LIMIT)


@pytest.fixture(scope="module")
def rust_binary():
    """Build Rust binary with netcdf support."""
    result = subprocess.run(
        ["cargo", "build", "--release", "--features", "netcdf", "--example", "process_pace"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        pytest.skip(f"Rust build failed: {result.stderr[:200]}")
    binary = REPO_ROOT / "target" / "release" / "examples" / "process_pace"
    assert binary.exists()
    return binary


@pytest.fixture(scope="module")
def rust_output(rust_binary):
    """Run Rust processing on cached PACE file."""
    if not PACE_FILE.exists():
        pytest.skip("No cached PACE file. Run: python tests/regression/run_pace_real.py")
    out_dir = CACHE_DIR / "rust_regression"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Clean previous output
    zarr_dir = out_dir / "PACE_OCI.20240701T175112.L1B.V3_corrected.zarr"
    if zarr_dir.exists():
        import shutil
        shutil.rmtree(zarr_dir)

    t0 = time.time()
    result = subprocess.run(
        [str(rust_binary), "--file", str(PACE_FILE),
         "--output", str(out_dir), "--limit", LIMIT_STR],
        capture_output=True, text=True, timeout=120,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0
    assert result.returncode == 0, f"Rust failed:\n{result.stdout}\n{result.stderr}"
    return {"stdout": result.stdout, "elapsed": elapsed, "dir": out_dir}


@pytest.fixture(scope="module")
def python_output():
    """Run Python ACOLITE processing on cached PACE file."""
    if not PACE_FILE.exists():
        pytest.skip("No cached PACE file")
    sys.path.insert(0, str(REPO_ROOT))
    try:
        import acolite as ac
        from netCDF4 import Dataset
    except ImportError:
        pytest.skip("acolite or netCDF4 not importable")

    out_dir = CACHE_DIR / "python_regression"
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    ofiles, _ = ac.pace.l1_convert(
        str(PACE_FILE), output=str(out_dir),
        settings={"limit": LIMIT, "verbosity": 0},
    )
    elapsed = time.time() - t0
    assert ofiles, "Python produced no output"

    # Extract stats
    stats = {}
    with Dataset(ofiles[0]) as nc:
        for v in sorted(nc.variables):
            if v.startswith("rhot_"):
                data = nc.variables[v][:]
                if hasattr(data, "mask"):
                    data = np.where(data.mask, np.nan, data.data)
                stats[v] = {
                    "mean": float(np.nanmean(data)),
                    "std": float(np.nanstd(data)),
                    "min": float(np.nanmin(data)),
                    "max": float(np.nanmax(data)),
                    "shape": list(data.shape),
                }
    return {"stats": stats, "elapsed": elapsed, "file": ofiles[0]}


def _read_zarr_v3_chunk(chunk_path, zarr_json_path):
    """Read a zarr v3 chunk, handling gzip or raw encoding."""
    meta = json.loads(Path(zarr_json_path).read_text())
    codecs = meta.get("codecs", [])
    has_gzip = any(c.get("name") == "gzip" for c in codecs)

    raw = Path(chunk_path).read_bytes()
    if has_gzip:
        import gzip
        raw = gzip.decompress(raw)
    return np.frombuffer(raw, dtype="<f4")


@pytest.fixture(scope="module")
def rust_stats(rust_output):
    """Extract per-band stats from Rust GeoZarr v3 output."""
    zarr_paths = list(rust_output["dir"].glob("*.zarr"))
    assert zarr_paths, "No zarr output found"
    zarr_dir = zarr_paths[0]

    # Read zarr v3 metadata directly
    data_meta = json.loads((zarr_dir / "data" / "zarr.json").read_text())
    shape = data_meta["shape"]  # [291, nrows, ncols]
    nbands, nrows, ncols = shape

    # Read wavelengths
    wavelengths = _read_zarr_v3_chunk(
        zarr_dir / "wavelengths" / "c" / "0",
        zarr_dir / "wavelengths" / "zarr.json",
    )

    # Read each band chunk
    stats = {}
    for i in range(nbands):
        chunk_path = zarr_dir / "data" / "c" / str(i) / "0" / "0"
        if not chunk_path.exists():
            continue
        band = _read_zarr_v3_chunk(
            chunk_path, zarr_dir / "data" / "zarr.json"
        ).reshape(nrows, ncols)
        valid = band[np.isfinite(band)]
        if len(valid) == 0:
            continue
        wl = float(wavelengths[i]) if i < len(wavelengths) else 0.0
        stats[f"band_{i}_wl{wl:.0f}"] = {
            "mean": float(np.mean(valid)),
            "std": float(np.std(valid)),
            "min": float(np.min(valid)),
            "max": float(np.max(valid)),
            "wavelength": wl,
            "shape": list(band.shape),
        }
    return stats


# ---------------------------------------------------------------------------
# Performance comparison tests
# ---------------------------------------------------------------------------

class TestPerformance:
    """Compare Rust vs Python processing performance."""

    def test_rust_faster_than_python(self, rust_output, python_output):
        """Rust total time must be faster than Python."""
        rust_t = rust_output["elapsed"]
        py_t = python_output["elapsed"]
        speedup = py_t / rust_t if rust_t > 0 else float("inf")
        print(f"\n  Python: {py_t:.2f}s")
        print(f"  Rust:   {rust_t:.2f}s")
        print(f"  Speedup: {speedup:.1f}×")
        assert rust_t < py_t, f"Rust ({rust_t:.2f}s) not faster than Python ({py_t:.2f}s)"

    def test_rust_ac_throughput(self, rust_output):
        """Rust AC throughput must exceed 10 Mpx/s."""
        stdout = rust_output["stdout"]
        # Parse "291 bands in 1.92ms (47.8 Mpx/s)"
        for line in stdout.split("\n"):
            if "Mpx/s" in line:
                mpx = float(line.split("(")[1].split("Mpx")[0].strip())
                print(f"\n  Rust AC throughput: {mpx:.1f} Mpx/s")
                assert mpx > 10, f"Throughput {mpx:.1f} Mpx/s below 10 Mpx/s threshold"
                return
        pytest.fail("Could not parse Mpx/s from Rust output")

    def test_rust_load_time(self, rust_output):
        """Rust load time for ROI subset must be under 5s."""
        stdout = rust_output["stdout"]
        for line in stdout.split("\n"):
            if "Loaded in" in line:
                # Parse "Loaded in 217.57ms" or "Loaded in 1.23s"
                time_str = line.split("Loaded in")[1].strip()
                if "ms" in time_str:
                    load_ms = float(time_str.replace("ms", ""))
                elif "s" in time_str:
                    load_ms = float(time_str.replace("s", "")) * 1000
                else:
                    continue
                print(f"\n  Rust load time: {load_ms:.0f}ms")
                assert load_ms < 5000, f"Load time {load_ms:.0f}ms exceeds 5s"
                return
        pytest.fail("Could not parse load time from Rust output")

    def test_rust_write_time(self, rust_output):
        """GeoZarr write time for 291 bands must be under 5s."""
        stdout = rust_output["stdout"]
        for line in stdout.split("\n"):
            if "Written in" in line:
                time_str = line.split("Written in")[1].strip()
                if "ms" in time_str:
                    write_ms = float(time_str.replace("ms", ""))
                elif "s" in time_str:
                    write_ms = float(time_str.replace("s", "")) * 1000
                else:
                    continue
                print(f"\n  Rust write time: {write_ms:.0f}ms")
                assert write_ms < 5000, f"Write time {write_ms:.0f}ms exceeds 5s"
                return
        pytest.fail("Could not parse write time from Rust output")


# ---------------------------------------------------------------------------
# Output comparison tests
# ---------------------------------------------------------------------------

class TestOutputComparison:
    """Compare Rust and Python output values."""

    def test_band_count_match(self, rust_stats, python_output):
        """Both must produce the same number of bands."""
        py_count = len(python_output["stats"])
        rust_count = len(rust_stats)
        print(f"\n  Python bands: {py_count}")
        print(f"  Rust bands:   {rust_count}")
        assert rust_count == py_count, f"Band count mismatch: Rust={rust_count} vs Python={py_count}"

    def test_wavelength_coverage(self, rust_stats, python_output):
        """Rust must cover the same wavelength range as Python."""
        py_bands = python_output["stats"]
        # Extract wavelengths from Python band names (e.g., rhot_blue_443 → 443)
        py_wls = sorted(set(
            int(k.rsplit("_", 1)[-1]) for k in py_bands.keys()
        ))
        rust_wls = sorted(set(
            int(round(v["wavelength"])) for v in rust_stats.values()
        ))
        print(f"\n  Python wavelength range: {py_wls[0]}–{py_wls[-1]} nm ({len(py_wls)} bands)")
        print(f"  Rust wavelength range:   {rust_wls[0]}–{rust_wls[-1]} nm ({len(rust_wls)} bands)")
        assert abs(py_wls[0] - rust_wls[0]) <= 2, "Start wavelength mismatch"
        assert abs(py_wls[-1] - rust_wls[-1]) <= 2, "End wavelength mismatch"

    def test_rhot_means_correlated(self, rust_stats, python_output):
        """Per-band mean reflectance must be correlated (r > 0.95).

        Note: Rust applies gas+Rayleigh+aerosol correction while Python L1R
        is TOA reflectance, so absolute values differ. But spectral shape
        (relative ordering) must be preserved.
        """
        py_bands = python_output["stats"]
        # Build wavelength→mean maps
        py_wl_mean = {}
        for k, v in py_bands.items():
            wl = int(k.rsplit("_", 1)[-1])
            py_wl_mean[wl] = v["mean"]

        rust_wl_mean = {}
        for v in rust_stats.values():
            wl = int(round(v["wavelength"]))
            rust_wl_mean[wl] = v["mean"]

        # Match by wavelength
        common_wls = sorted(set(py_wl_mean.keys()) & set(rust_wl_mean.keys()))
        assert len(common_wls) > 200, f"Only {len(common_wls)} common wavelengths"

        py_vals = np.array([py_wl_mean[w] for w in common_wls])
        rust_vals = np.array([rust_wl_mean[w] for w in common_wls])

        corr = np.corrcoef(py_vals, rust_vals)[0, 1]
        print(f"\n  Common wavelengths: {len(common_wls)}")
        print(f"  Correlation (Python TOA vs Rust corrected): {corr:.4f}")
        assert corr > 0.95, f"Correlation {corr:.4f} below 0.95 threshold"

    def test_rust_corrected_lower_than_toa(self, rust_stats, python_output):
        """Rust corrected reflectance should generally be ≤ Python TOA.

        Atmospheric correction removes Rayleigh/aerosol, so ρs ≤ ρt for most bands.
        Allow some bands to be slightly higher due to correction artifacts.
        """
        py_bands = python_output["stats"]
        py_wl_mean = {int(k.rsplit("_", 1)[-1]): v["mean"] for k, v in py_bands.items()}
        rust_wl_mean = {int(round(v["wavelength"])): v["mean"] for v in rust_stats.values()}

        common = sorted(set(py_wl_mean.keys()) & set(rust_wl_mean.keys()))
        lower_count = sum(1 for w in common if rust_wl_mean[w] <= py_wl_mean[w] + 0.01)
        pct = lower_count / len(common) * 100
        print(f"\n  Bands where Rust ρs ≤ Python ρt (+0.01): {lower_count}/{len(common)} ({pct:.0f}%)")
        assert pct > 70, f"Only {pct:.0f}% of bands have Rust ≤ Python"

    def test_spectral_shape_preserved(self, rust_stats):
        """Rust output must preserve spectral features (O2-A dip, H2O absorption)."""
        wl_mean = {int(round(v["wavelength"])): v["mean"] for v in rust_stats.values()}

        # O2-A absorption dip at ~761 nm
        if all(w in wl_mean for w in [754, 761, 771]):
            assert wl_mean[761] < wl_mean[754], "No O2-A dip in Rust output"
            assert wl_mean[761] < wl_mean[771], "No O2-A dip in Rust output"
            print(f"\n  O2-A dip: ρs(754)={wl_mean[754]:.4f} > ρs(761)={wl_mean[761]:.4f} < ρs(771)={wl_mean[771]:.4f} ✓")

        # Blue > SWIR (Rayleigh)
        blue_mean = np.mean([v for w, v in wl_mean.items() if 400 <= w <= 500])
        swir_mean = np.mean([v for w, v in wl_mean.items() if w >= 1200])
        print(f"  Blue mean: {blue_mean:.4f}, SWIR mean: {swir_mean:.4f}")
        assert blue_mean > swir_mean, "Blue not > SWIR in Rust output"


# ---------------------------------------------------------------------------
# GeoZarr structure tests
# ---------------------------------------------------------------------------

class TestGeoZarrOutput:
    """Validate Rust GeoZarr output structure."""

    def test_zarr_exists(self, rust_output):
        zarr_dirs = list(rust_output["dir"].glob("*.zarr"))
        assert len(zarr_dirs) == 1, f"Expected 1 zarr dir, found {len(zarr_dirs)}"

    def test_zarr_metadata(self, rust_output):
        zarr_dir = list(rust_output["dir"].glob("*.zarr"))[0]
        meta_path = zarr_dir / "zarr.json"
        assert meta_path.exists(), "Missing zarr.json"
        meta = json.loads(meta_path.read_text())
        attrs = meta.get("attributes", {})
        assert attrs.get("Conventions") == "CF-1.8"
        assert attrs.get("sensor") == "PACE_OCI"
        assert "proj:epsg" in attrs

    def test_zarr_data_shape(self, rust_output):
        zarr_dir = list(rust_output["dir"].glob("*.zarr"))[0]
        meta = json.loads((zarr_dir / "data" / "zarr.json").read_text())
        shape = meta["shape"]
        assert len(shape) == 3, f"Expected 3D, got {len(shape)}D"
        assert shape[0] == 291, f"Expected 291 bands, got {shape[0]}"
        print(f"\n  GeoZarr shape: {shape}")

    def test_zarr_wavelengths(self, rust_output):
        zarr_dir = list(rust_output["dir"].glob("*.zarr"))[0]
        wl = _read_zarr_v3_chunk(
            zarr_dir / "wavelengths" / "c" / "0",
            zarr_dir / "wavelengths" / "zarr.json",
        )
        assert len(wl) == 291
        assert wl[0] >= 300 and wl[-1] <= 2300
        # Blue and red detectors are internally sorted; SWIR has near-duplicate
        # bands (1250/1249, 1620/1618) so only check overall range
        blue_wl = wl[:119]
        red_wl = wl[119:119+163]
        swir_wl = wl[119+163:]
        assert np.all(np.diff(blue_wl) > 0), "Blue wavelengths not sorted"
        assert np.all(np.diff(red_wl) > 0), "Red wavelengths not sorted"
        assert swir_wl[0] >= 900 and swir_wl[-1] <= 2300, "SWIR range wrong"


# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

class TestRegressionSummary:
    """Generate summary report."""

    def test_print_summary(self, rust_output, python_output, rust_stats):
        py_t = python_output["elapsed"]
        rust_t = rust_output["elapsed"]
        speedup = py_t / rust_t if rust_t > 0 else float("inf")

        print(f"\n{'='*60}")
        print("PACE OCI REGRESSION SUMMARY")
        print(f"{'='*60}")
        print(f"  Granule: PACE_OCI.20240701T175112.L1B.V3.nc")
        print(f"  ROI: {LIMIT}")
        print(f"  Bands: {len(rust_stats)}")
        print(f"")
        print(f"  Python ACOLITE (L1R):  {py_t:.2f}s")
        print(f"  Rust acolite-rs (AC):  {rust_t:.2f}s")
        print(f"  Speedup:               {speedup:.1f}×")
        print(f"")

        # Parse Rust timing breakdown
        for line in rust_output["stdout"].split("\n"):
            if any(k in line for k in ["Loaded in", "bands in", "Written in"]):
                print(f"  {line.strip()}")
        print(f"{'='*60}")
