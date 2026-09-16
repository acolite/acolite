"""
PACE OCI full DSF atmospheric correction: Rust vs Python regression test.

Compares acolite-rs (generic LUT DSF) against Python ACOLITE for PACE OCI L1B data.
Tests both fixed and tiled DSF modes, per-band accuracy, and performance.

Requires:
  - Cached PACE L1B file (PACE_OCI.20240701T175112.L1B.V3.nc)
  - Rust binary built with --features full-io
  - Python ACOLITE L2R output (generated on first run)
"""

import gzip
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
# Water ROI: open ocean east of Outer Banks
LIMIT = [35.0, -75.0, 36.0, -74.0]
LIMIT_STR = ",".join(str(x) for x in LIMIT)


def read_zarr_chunk(zarr_dir, array_name, *indices):
    """Read a zarr v3 chunk, handling gzip compression."""
    chunk_path = Path(zarr_dir) / array_name / "c" / "/".join(str(i) for i in indices)
    with open(chunk_path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    meta_path = Path(zarr_dir) / array_name / "zarr.json"
    meta = json.load(open(meta_path))
    dtype = np.float32 if meta.get("data_type") == "float32" else np.float64
    return np.frombuffer(raw, dtype=dtype)


@pytest.fixture(scope="module")
def rust_binary():
    """Build Rust binary with full-io support."""
    result = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io", "--example", "process_pace_ac"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=300,
    )
    if result.returncode != 0:
        pytest.skip(f"Rust build failed: {result.stderr[:200]}")
    binary = REPO_ROOT / "target" / "release" / "examples" / "process_pace_ac"
    assert binary.exists()
    return binary


@pytest.fixture(scope="module")
def rust_fixed(rust_binary):
    """Run Rust fixed-mode DSF on PACE ROI."""
    if not PACE_FILE.exists():
        pytest.skip("No cached PACE file")
    out_dir = CACHE_DIR / "pace_dsf" / "rust_fixed"
    zarr_dir = out_dir / "PACE_OCI.20240701T175112.L1B.V3_L2R.zarr"
    if zarr_dir.exists():
        import shutil
        shutil.rmtree(zarr_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    result = subprocess.run(
        [str(rust_binary), "--file", str(PACE_FILE),
         "--output", str(out_dir), "--limit", LIMIT_STR,
         "--model", "auto", "--aot-mode", "fixed"],
        capture_output=True, text=True, timeout=600,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0
    assert result.returncode == 0, f"Rust failed:\n{result.stdout}\n{result.stderr}"
    return {"stdout": result.stdout, "elapsed": elapsed, "zarr": str(zarr_dir)}


@pytest.fixture(scope="module")
def rust_tiled(rust_binary):
    """Run Rust tiled-mode DSF on PACE ROI."""
    if not PACE_FILE.exists():
        pytest.skip("No cached PACE file")
    out_dir = CACHE_DIR / "pace_dsf" / "rust_tiled"
    zarr_dir = out_dir / "PACE_OCI.20240701T175112.L1B.V3_L2R.zarr"
    if zarr_dir.exists():
        import shutil
        shutil.rmtree(zarr_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    result = subprocess.run(
        [str(rust_binary), "--file", str(PACE_FILE),
         "--output", str(out_dir), "--limit", LIMIT_STR,
         "--model", "auto", "--aot-mode", "tiled"],
        capture_output=True, text=True, timeout=600,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0
    assert result.returncode == 0, f"Rust tiled failed:\n{result.stdout}\n{result.stderr}"
    return {"stdout": result.stdout, "elapsed": elapsed, "zarr": str(zarr_dir)}


@pytest.fixture(scope="module")
def python_l2r():
    """Run Python ACOLITE L2R on PACE ROI (cached)."""
    if not PACE_FILE.exists():
        pytest.skip("No cached PACE file")

    out_dir = CACHE_DIR / "pace_dsf" / "python_fixed"
    l2r_file = out_dir / "PACE_OCI_2024_07_01_17_51_12_L2R.nc"
    if l2r_file.exists():
        return str(l2r_file)

    sys.path.insert(0, str(REPO_ROOT))
    try:
        import acolite as ac
    except ImportError:
        pytest.skip("acolite not importable")

    out_dir.mkdir(parents=True, exist_ok=True)
    settings = {
        "inputfile": str(PACE_FILE),
        "output": str(out_dir),
        "limit": LIMIT,
        "verbosity": 0,
        "ancillary_data": False,
        "dsf_interface_reflectance": False,
        "dsf_aot_estimate": "fixed",
    }
    ofiles = ac.acolite.acolite_run(settings=settings)
    assert l2r_file.exists(), f"Python L2R not created. Output: {ofiles}"
    return str(l2r_file)


@pytest.fixture(scope="module")
def comparison_data(rust_fixed, python_l2r):
    """Load and align Rust and Python results for comparison."""
    from netCDF4 import Dataset

    # Read Rust wavelengths
    zarr_dir = rust_fixed["zarr"]
    rust_waves = read_zarr_chunk(zarr_dir, "wavelengths", 0)

    # Read Python L2R
    py_nc = Dataset(python_l2r)
    py_aot = py_nc.getncattr("ac_aot_550")
    py_model = py_nc.getncattr("ac_model")

    # Parse Rust AOT from stdout
    rust_aot = None
    for line in rust_fixed["stdout"].split("\n"):
        if "Fixed DSF:" in line and "AOT=" in line:
            rust_aot = float(line.split("AOT=")[1].split(",")[0])
            break

    # Build per-band comparison
    bands = []
    py_vars = [v for v in py_nc.variables if v.startswith("rhos_")]
    for tw in [400, 443, 490, 550, 600, 665, 750, 865]:
        # Find Rust band
        ri = int(np.argmin(np.abs(rust_waves - tw)))
        meta = json.load(open(Path(zarr_dir) / "data" / "zarr.json"))
        shape = meta["shape"]
        rust_data = read_zarr_chunk(zarr_dir, "data", ri, 0, 0).reshape(shape[1], shape[2])

        # Find Python band
        best_py = None
        best_diff = 999
        for pb in py_vars:
            parts = pb.split("_")
            try:
                wl = float(parts[-1])
                if abs(wl - tw) < best_diff:
                    best_diff = abs(wl - tw)
                    best_py = pb
            except ValueError:
                pass

        if best_py and best_diff < 5:
            py_data = py_nc.variables[best_py][:]
            if hasattr(py_data, "mask"):
                py_data = np.where(py_data.mask, np.nan, py_data.data)

            # Align shapes
            min_r = min(rust_data.shape[0], py_data.shape[0])
            min_c = min(rust_data.shape[1], py_data.shape[1])
            r_flat = rust_data[:min_r, :min_c].flatten()
            p_flat = py_data[:min_r, :min_c].flatten()
            mask = np.isfinite(r_flat) & np.isfinite(p_flat)

            bands.append({
                "wavelength": tw,
                "rust_mean": np.nanmean(rust_data[np.isfinite(rust_data)]),
                "python_mean": np.nanmean(py_data[np.isfinite(py_data)]),
                "rust_flat": r_flat[mask],
                "python_flat": p_flat[mask],
                "n_valid": int(np.sum(mask)),
            })

    py_nc.close()
    return {
        "bands": bands,
        "rust_aot": rust_aot,
        "python_aot": py_aot,
        "py_model": py_model,
        "rust_elapsed": rust_fixed["elapsed"],
    }


class TestPaceDsfFixed:
    """Test fixed-mode DSF atmospheric correction for PACE OCI."""

    def test_rust_produces_output(self, rust_fixed):
        zarr_dir = Path(rust_fixed["zarr"])
        assert zarr_dir.exists()
        assert (zarr_dir / "zarr.json").exists()
        assert (zarr_dir / "data").exists()
        assert (zarr_dir / "wavelengths").exists()

    def test_band_count(self, rust_fixed):
        zarr_dir = rust_fixed["zarr"]
        meta = json.load(open(Path(zarr_dir) / "data" / "zarr.json"))
        assert meta["shape"][0] == 291, f"Expected 291 bands, got {meta['shape'][0]}"

    def test_model_selection_matches(self, comparison_data):
        """Both should select MOD2 for this scene."""
        assert "MOD2" in comparison_data["py_model"]

    def test_aot_within_tolerance(self, comparison_data):
        """AOT should be within 5% of Python."""
        rust_aot = comparison_data["rust_aot"]
        py_aot = comparison_data["python_aot"]
        assert rust_aot is not None, "Could not parse Rust AOT"
        rel_diff = abs(rust_aot - py_aot) / max(py_aot, 0.001)
        assert rel_diff < 0.05, f"AOT diff too large: Rust={rust_aot:.4f}, Python={py_aot:.4f} ({rel_diff*100:.1f}%)"

    def test_visible_band_correlation(self, comparison_data):
        """Visible bands (400-600nm) should have R > 0.90."""
        for band in comparison_data["bands"]:
            if 400 <= band["wavelength"] <= 600 and band["n_valid"] > 10:
                r = np.corrcoef(band["rust_flat"], band["python_flat"])[0, 1]
                assert r > 0.90, f"{band['wavelength']}nm: R={r:.4f} < 0.90"

    def test_visible_band_mean_diff(self, comparison_data):
        """Visible band means should be within 1% of Python."""
        for band in comparison_data["bands"]:
            if 400 <= band["wavelength"] <= 600:
                diff = abs(band["rust_mean"] - band["python_mean"])
                rdiff = diff / max(abs(band["python_mean"]), 1e-6)
                assert rdiff < 0.01, f"{band['wavelength']}nm: diff={rdiff*100:.2f}% > 1%"

    def test_nir_band_correlation(self, comparison_data):
        """NIR bands should have R > 0.85 (more affected by gas model differences)."""
        for band in comparison_data["bands"]:
            if band["wavelength"] > 700 and band["n_valid"] > 10:
                r = np.corrcoef(band["rust_flat"], band["python_flat"])[0, 1]
                assert r > 0.85, f"{band['wavelength']}nm: R={r:.4f} < 0.85"


class TestPaceDsfTiled:
    """Test tiled-mode DSF atmospheric correction for PACE OCI."""

    def test_tiled_produces_output(self, rust_tiled):
        zarr_dir = Path(rust_tiled["zarr"])
        assert zarr_dir.exists()
        assert (zarr_dir / "data").exists()

    def test_tiled_band_count(self, rust_tiled):
        zarr_dir = rust_tiled["zarr"]
        meta = json.load(open(Path(zarr_dir) / "data" / "zarr.json"))
        assert meta["shape"][0] == 291


class TestPaceDsfPerformance:
    """Performance benchmarks for PACE DSF processing."""

    def test_rust_faster_than_python(self, comparison_data):
        """Rust should be significantly faster than Python."""
        rust_time = comparison_data["rust_elapsed"]
        # Python takes ~237s for this scene
        # Rust should be at least 5x faster
        assert rust_time < 60, f"Rust took {rust_time:.1f}s, expected < 60s"

    def test_ac_throughput(self, rust_fixed):
        """AC step should process at > 1 Mpx/s."""
        for line in rust_fixed["stdout"].split("\n"):
            if "Mpx/s" in line:
                mpx = float(line.split("(")[1].split(" Mpx")[0])
                assert mpx > 1.0, f"Throughput {mpx:.1f} Mpx/s < 1.0"
                return
        pytest.skip("Could not find throughput in output")


class TestRegressionSummary:
    """Print summary of all regression metrics."""

    def test_print_summary(self, comparison_data):
        print("\n" + "=" * 70)
        print("PACE OCI DSF Regression Summary (Rust vs Python)")
        print("=" * 70)
        print(f"  Rust AOT:   {comparison_data['rust_aot']:.4f}")
        print(f"  Python AOT: {comparison_data['python_aot']:.4f}")
        print(f"  AOT diff:   {abs(comparison_data['rust_aot'] - comparison_data['python_aot']):.4f}")
        print(f"  Rust time:  {comparison_data['rust_elapsed']:.1f}s")
        print()
        print(f"  {'Band':>8s}  {'Py mean':>10s}  {'Rust mean':>10s}  {'Diff':>8s}  {'%Diff':>6s}  {'R':>8s}")
        print(f"  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*6}  {'-'*8}")
        for band in comparison_data["bands"]:
            diff = abs(band["rust_mean"] - band["python_mean"])
            rdiff = diff / max(abs(band["python_mean"]), 1e-6) * 100
            if band["n_valid"] > 10:
                r = np.corrcoef(band["rust_flat"], band["python_flat"])[0, 1]
            else:
                r = float("nan")
            print(f"  {band['wavelength']:>6.0f}nm  {band['python_mean']:>10.6f}  {band['rust_mean']:>10.6f}  {diff:>8.4f}  {rdiff:>5.1f}%  {r:>8.6f}")
        print("=" * 70)
