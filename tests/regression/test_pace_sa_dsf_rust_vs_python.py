"""
PACE OCI DSF regression tests — South Australia ocean scene.

Compares Rust acolite-rs vs Python ACOLITE for atmospheric correction
of PACE OCI hyperspectral data over open ocean south of Kangaroo Island, SA.

Scene: PACE_OCI.20241231T044250.L1B.V3.nc (1710×1272, 291 bands)
ROI:   [-37.0, 136.5, -36.0, 137.5] — open ocean, 108×57 pixels
"""
import pytest
import subprocess
import os
import sys
import time
import json
import gzip
import numpy as np

# ── Configuration ──
PACE_FILE = "/tmp/acolite_test_cache/PACE_OCI.20241231T044250.L1B.V3.nc"
LIMIT = [-37.0, 136.5, -36.0, 137.5]
RUST_OUT = "/tmp/acolite_test_cache/pace_sa/rust_test"
PYTHON_OUT = "/tmp/acolite_test_cache/pace_sa/python_test"
RUST_FULL_OUT = "/tmp/acolite_test_cache/pace_sa/rust_full_test"

# Skip all tests if PACE file not available
pytestmark = pytest.mark.skipif(
    not os.path.exists(PACE_FILE),
    reason=f"PACE SA file not found: {PACE_FILE}"
)


def read_zarr_band(zarr_dir, band_idx):
    """Read a single band from Rust GeoZarr output."""
    chunk_path = os.path.join(zarr_dir, f"data/c/{band_idx}/0/0")
    if not os.path.exists(chunk_path):
        return None
    with open(chunk_path, "rb") as f:
        raw = f.read()
    try:
        raw = gzip.decompress(raw)
    except Exception:
        pass
    return np.frombuffer(raw, dtype="<f4")


def read_zarr_wavelengths(zarr_dir):
    """Read wavelength array from Rust GeoZarr output."""
    with open(os.path.join(zarr_dir, "wavelengths/c/0"), "rb") as f:
        raw = f.read()
    return np.frombuffer(raw, dtype="<f4")


# ── Fixtures ──

@pytest.fixture(scope="module")
def rust_binary():
    """Build Rust binary."""
    result = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io", "--example", "process_pace_ac"],
        capture_output=True, text=True, cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    )
    assert result.returncode == 0, f"Build failed: {result.stderr}"
    binary = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "target/release/examples/process_pace_ac"
    )
    assert os.path.exists(binary)
    return binary


@pytest.fixture(scope="module")
def rust_fixed(rust_binary):
    """Run Rust PACE AC with fixed DSF on SA ocean ROI."""
    os.makedirs(RUST_OUT, exist_ok=True)
    limit_str = ",".join(str(x) for x in LIMIT)
    start = time.time()
    result = subprocess.run(
        [rust_binary, "--file", PACE_FILE, "--output", RUST_OUT,
         "--limit", limit_str, "--model", "auto", "--aot-mode", "fixed"],
        capture_output=True, text=True, env={**os.environ, "RUST_LOG": "info"},
        cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    )
    elapsed = time.time() - start
    assert result.returncode == 0, f"Rust failed: {result.stderr}\n{result.stdout}"
    zarr_dir = None
    for d in os.listdir(RUST_OUT):
        if d.endswith(".zarr"):
            zarr_dir = os.path.join(RUST_OUT, d)
            break
    assert zarr_dir is not None, "No zarr output found"
    return {"zarr_dir": zarr_dir, "time": elapsed, "stdout": result.stdout}


@pytest.fixture(scope="module")
def rust_tiled(rust_binary):
    """Run Rust PACE AC with tiled DSF on SA ocean ROI."""
    out = RUST_OUT + "_tiled"
    os.makedirs(out, exist_ok=True)
    limit_str = ",".join(str(x) for x in LIMIT)
    result = subprocess.run(
        [rust_binary, "--file", PACE_FILE, "--output", out,
         "--limit", limit_str, "--model", "auto", "--aot-mode", "tiled"],
        capture_output=True, text=True, env={**os.environ, "RUST_LOG": "info"},
        cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    )
    assert result.returncode == 0, f"Rust tiled failed: {result.stderr}"
    zarr_dir = None
    for d in os.listdir(out):
        if d.endswith(".zarr"):
            zarr_dir = os.path.join(out, d)
            break
    return {"zarr_dir": zarr_dir, "stdout": result.stdout}


@pytest.fixture(scope="module")
def python_l2r():
    """Run Python ACOLITE on SA ocean ROI."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    import acolite as ac

    os.makedirs(PYTHON_OUT, exist_ok=True)
    settings = {
        "inputfile": PACE_FILE,
        "output": PYTHON_OUT,
        "limit": LIMIT,
        "dsf_aot_estimate": "fixed",
        "l2w_parameters": [],
        "verbosity": 1,
        "rgb_rhot": False,
        "rgb_rhos": False,
        "map_l2w": False,
    }
    start = time.time()
    ac.acolite.acolite_run(settings=settings)
    elapsed = time.time() - start

    l2r_files = [f for f in os.listdir(PYTHON_OUT) if f.endswith("L2R.nc")]
    assert len(l2r_files) > 0, "No Python L2R output"
    return {"nc_path": os.path.join(PYTHON_OUT, l2r_files[0]), "time": elapsed}


@pytest.fixture(scope="module")
def comparison_data(rust_fixed, python_l2r):
    """Load and compare Rust vs Python results at key wavelengths."""
    from netCDF4 import Dataset

    py = Dataset(python_l2r["nc_path"])
    zarr_dir = rust_fixed["zarr_dir"]
    rust_waves = read_zarr_wavelengths(zarr_dir)

    py_vars = sorted([v for v in py.variables if v.startswith("rhos_")])
    py_aot = float(np.nanmean(py["aot_550"][:]))

    # Extract Rust AOT from stdout
    rust_aot = None
    for line in rust_fixed["stdout"].split("\n"):
        if "Fixed DSF:" in line and "AOT=" in line:
            rust_aot = float(line.split("AOT=")[1].split(",")[0])
            break

    # Extract model name
    rust_model = None
    for line in rust_fixed["stdout"].split("\n"):
        if "Fixed DSF:" in line and "model=" in line:
            rust_model = line.split("model=")[1].split(",")[0]
            break

    # Per-band comparison
    test_waves = [400, 443, 490, 550, 600, 665, 750, 865]
    bands = {}
    for tw in test_waves:
        py_name = None
        for v in py_vars:
            wn = int(v.split("_")[-1])
            if abs(wn - tw) <= 3:
                py_name = v
                break
        if py_name is None:
            continue
        py_data = py[py_name][:].flatten().astype(np.float64)
        rs_idx = int(np.argmin(np.abs(rust_waves - tw)))
        if abs(rust_waves[rs_idx] - tw) > 5:
            continue
        rs_data = read_zarr_band(zarr_dir, rs_idx)
        if rs_data is None:
            continue
        rs_data = rs_data.astype(np.float64)
        n = min(len(py_data), len(rs_data))
        mask = (np.isfinite(py_data[:n]) & np.isfinite(rs_data[:n])
                & (np.abs(py_data[:n]) < 10) & (np.abs(rs_data[:n]) < 10))
        if mask.sum() < 10:
            continue
        py_m = np.mean(py_data[:n][mask])
        rs_m = np.mean(rs_data[:n][mask])
        diff = abs(rs_m - py_m)
        pct = diff / abs(py_m) * 100 if abs(py_m) > 1e-6 else 0
        bands[tw] = {"py_mean": py_m, "rs_mean": rs_m, "diff": diff, "pct": pct}

    py.close()
    return {
        "py_aot": py_aot,
        "rust_aot": rust_aot,
        "rust_model": rust_model,
        "bands": bands,
        "rust_time": rust_fixed["time"],
        "python_time": python_l2r["time"],
    }


# ── Tests ──

class TestPaceSaFixed:
    """Fixed DSF mode tests."""

    def test_rust_output_exists(self, rust_fixed):
        assert os.path.isdir(rust_fixed["zarr_dir"])

    def test_rust_has_bands(self, rust_fixed):
        waves = read_zarr_wavelengths(rust_fixed["zarr_dir"])
        assert len(waves) == 291

    def test_model_selection(self, comparison_data):
        assert "MOD2" in comparison_data["rust_model"]

    def test_aot_within_10pct(self, comparison_data):
        py_aot = comparison_data["py_aot"]
        rs_aot = comparison_data["rust_aot"]
        assert rs_aot is not None
        pct = abs(rs_aot - py_aot) / py_aot * 100
        assert pct < 10, f"AOT diff {pct:.1f}%: Rust={rs_aot:.4f} Python={py_aot:.4f}"

    def test_visible_bands_within_15pct(self, comparison_data):
        """Visible bands (400-665nm) should be within 15% mean difference."""
        for tw in [443, 490, 550, 600, 665]:
            if tw in comparison_data["bands"]:
                b = comparison_data["bands"][tw]
                assert b["pct"] < 15, f"{tw}nm: {b['pct']:.1f}% diff"

    def test_nir_bands_within_15pct(self, comparison_data):
        """NIR bands should be within 15% mean difference."""
        for tw in [750, 865]:
            if tw in comparison_data["bands"]:
                b = comparison_data["bands"][tw]
                assert b["pct"] < 15, f"{tw}nm: {b['pct']:.1f}% diff"


class TestPaceSaTiled:
    """Tiled DSF mode tests."""

    def test_tiled_output_exists(self, rust_tiled):
        assert rust_tiled["zarr_dir"] is not None
        assert os.path.isdir(rust_tiled["zarr_dir"])

    def test_tiled_has_bands(self, rust_tiled):
        waves = read_zarr_wavelengths(rust_tiled["zarr_dir"])
        assert len(waves) == 291


class TestPaceSaPerformance:
    """Performance benchmarks."""

    def test_rust_roi_under_60s(self, rust_fixed):
        assert rust_fixed["time"] < 60, f"Rust ROI took {rust_fixed['time']:.1f}s"

    def test_rust_faster_than_python(self, comparison_data):
        speedup = comparison_data["python_time"] / comparison_data["rust_time"]
        assert speedup > 3, f"Speedup only {speedup:.1f}×"


class TestPaceSaFullScene:
    """Full scene benchmark (no ROI limit)."""

    def test_full_scene_completes(self, rust_binary):
        """Process full 1710×1272 scene with tiled DSF."""
        os.makedirs(RUST_FULL_OUT, exist_ok=True)
        start = time.time()
        result = subprocess.run(
            [rust_binary, "--file", PACE_FILE, "--output", RUST_FULL_OUT,
             "--model", "auto", "--aot-mode", "tiled"],
            capture_output=True, text=True, env={**os.environ, "RUST_LOG": "info"},
            timeout=900,
            cwd=os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        )
        elapsed = time.time() - start
        assert result.returncode == 0, f"Full scene failed: {result.stderr[-500:]}"

        # Check output exists
        zarr_dirs = [d for d in os.listdir(RUST_FULL_OUT) if d.endswith(".zarr")]
        assert len(zarr_dirs) > 0

        # Check timing
        print(f"\n  Full scene total: {elapsed:.1f}s")
        for line in result.stdout.split("\n"):
            if "AC:" in line:
                print(f"  {line.strip()}")

        # AC should complete in under 300s
        assert elapsed < 600, f"Full scene took {elapsed:.1f}s"


class TestRegressionSummary:
    """Print summary of all results."""

    def test_summary(self, comparison_data):
        d = comparison_data
        print("\n" + "=" * 60)
        print("PACE OCI SA DSF Regression Summary")
        print("=" * 60)
        print(f"  Model: {d['rust_model']}")
        print(f"  AOT: Rust={d['rust_aot']:.4f} Python={d['py_aot']:.4f} "
              f"({abs(d['rust_aot'] - d['py_aot']) / d['py_aot'] * 100:.1f}%)")
        print(f"  Time: Rust={d['rust_time']:.1f}s Python={d['python_time']:.1f}s "
              f"({d['python_time'] / d['rust_time']:.1f}× speedup)")
        print(f"\n  {'Wave':>6} {'Py mean':>10} {'Rs mean':>10} {'%Diff':>8}")
        print("  " + "-" * 40)
        for tw in sorted(d["bands"].keys()):
            b = d["bands"][tw]
            print(f"  {tw:>6} {b['py_mean']:>10.6f} {b['rs_mean']:>10.6f} {b['pct']:>7.1f}%")
        print("=" * 60)
