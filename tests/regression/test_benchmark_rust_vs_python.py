"""
Full-scene benchmark: Rust vs Python ACOLITE for Landsat 8/9.

Runs BOTH pipelines on the complete scene (no ROI subsetting) and measures:
  - Wall-clock time for each pipeline
  - Per-pixel accuracy (Pearson R, RMSE, bias)
  - Speedup factor

South Australia, Gulf St Vincent water scenes (7931×7891 pixels).
"""

import os, re, subprocess, sys, tempfile, time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent.parent
CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache")) / "landsat_sa"

BAND_MAP_L8 = {"B1": "443", "B2": "483", "B3": "561", "B4": "655", "B5": "865", "B6": "1609", "B7": "2201"}
BAND_MAP_L9 = {"B1": "443", "B2": "482", "B3": "561", "B4": "654", "B5": "865", "B6": "1608", "B7": "2201"}

SCENES = {
    "L8": {
        "id": "LC08_L1TP_098084_20240205_20240212_02_T1",
        "band_map": BAND_MAP_L8,
        "py_prefix": "L8_OLI_2024_02_05_00_39_50_098084",
    },
    "L9": {
        "id": "LC09_L1TP_098084_20240213_20240213_02_T1",
        "band_map": BAND_MAP_L9,
        "py_prefix": "L9_OLI_2024_02_13_00_40_00_098084",
    },
}


# ── Python ACOLITE full-scene processing ─────────────────────────────────

def _run_python_full(label):
    """Run Python ACOLITE on full scene (no limit/ROI)."""
    cfg = SCENES[label]
    scene_dir = CACHE / cfg["id"]
    out_dir = CACHE / f"py_full_{label}"

    if not scene_dir.exists():
        pytest.skip(f"No cached {label} scene")

    # Check for cached result
    l2r_nc = out_dir / f"{cfg['py_prefix']}_L2R.nc"
    if l2r_nc.exists():
        # Extract timing from log if available
        log = out_dir / "timing.txt"
        py_time = float(log.read_text().strip()) if log.exists() else None
        return {"nc": l2r_nc, "out_dir": out_dir, "time": py_time, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)

    # Write settings file — full scene, no limit
    # Use fixed AOT estimation (single AOT for whole scene) because
    # tiled estimation fails on full scenes with deep ocean tiles
    # (dark spectrum intercept produces values below LUT minimum romix).
    settings = out_dir / "settings.txt"
    settings.write_text(
        f"inputfile={scene_dir}\n"
        f"output={out_dir}\n"
        "dsf_aot_estimate=fixed\n"
        "dsf_spectrum_option=intercept\n"
        "dsf_intercept_pixels=200\n"
        "dsf_aot_compute=min\n"
        "dsf_nbands=2\n"
        "dsf_nbands_fit=2\n"
        "dsf_model_selection=min_drmsd\n"
        "dsf_interface_reflectance=False\n"
        "l2w_parameters=\n"
        "rgb_rhot=False\n"
        "rgb_rhos=False\n"
        "map_l2w=False\n"
    )

    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0,'{REPO}'); "
         f"import acolite; acolite.acolite.acolite_run(settings='{settings}')"],
        capture_output=True, text=True, timeout=600,
        cwd=str(REPO),
    )
    elapsed = time.time() - t0

    if r.returncode != 0:
        pytest.fail(f"Python {label} failed ({elapsed:.1f}s):\n{r.stderr[:1000]}")

    # Save timing
    (out_dir / "timing.txt").write_text(f"{elapsed:.2f}")

    l2r_nc = list(out_dir.glob("*_L2R.nc"))
    if not l2r_nc:
        pytest.fail(f"No L2R NetCDF produced for {label}:\n{r.stdout[-500:]}")

    return {"nc": l2r_nc[0], "out_dir": out_dir, "time": elapsed, "cached": False}


# ── Rust full-scene processing ───────────────────────────────────────────

def _run_rust_full(binary, label):
    """Run Rust on full scene."""
    cfg = SCENES[label]
    scene_dir = CACHE / cfg["id"]
    out_dir = CACHE / f"rust_full_{label}"

    if not scene_dir.exists():
        pytest.skip(f"No cached {label} scene")

    # Check for cached result
    tifs = sorted(out_dir.glob(f"{cfg['id']}_corrected_B*.tif")) if out_dir.exists() else []
    log = out_dir / "timing.txt"
    if len(tifs) == 7 and log.exists():
        rust_time = float(log.read_text().strip())
        return {"tifs": tifs, "out_dir": out_dir, "time": rust_time, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(scene_dir), "--output", str(out_dir)],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0

    if r.returncode != 0:
        pytest.fail(f"Rust {label} failed ({elapsed:.1f}s):\n{r.stderr[:500]}")

    # Parse sub-timings from stdout
    timings = {"total": elapsed}
    for line in r.stdout.splitlines():
        for key in ("Load", "AC", "Write"):
            m = re.search(rf"{key}:\s*([\d.]+)s", line)
            if m:
                timings[key.lower()] = float(m.group(1))

    (out_dir / "timing.txt").write_text(f"{elapsed:.2f}")

    tifs = sorted(out_dir.glob(f"{cfg['id']}_corrected_B*.tif"))
    return {"tifs": tifs, "out_dir": out_dir, "time": elapsed, "timings": timings, "cached": False}


# ── Pixel comparison helpers ─────────────────────────────────────────────

def _read_python_band(nc_path, wave_name):
    """Read rhos band from Python L2R NetCDF."""
    import netCDF4
    ds = netCDF4.Dataset(str(nc_path))
    var_name = f"rhos_{wave_name}"
    if var_name not in ds.variables:
        ds.close()
        return None
    data = ds.variables[var_name][:].astype(np.float64)
    ds.close()
    return data


def _read_rust_band(tif_path):
    """Read band from Rust COG GeoTIFF."""
    from osgeo import gdal
    ds = gdal.Open(str(tif_path))
    data = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    gt = ds.GetGeoTransform()
    ds = None
    return data, gt


def _align_and_compare(rust_data, rust_gt, py_data, py_shape):
    """Align Rust full-scene output to Python full-scene output and compute metrics."""
    # Both are full scene — should be same dimensions
    # Crop to common size if slightly different
    h = min(rust_data.shape[0], py_data.shape[0])
    w = min(rust_data.shape[1], py_data.shape[1])
    r = rust_data[:h, :w]
    p = py_data[:h, :w]

    mask = np.isfinite(r) & np.isfinite(p) & (np.abs(p) < 1e10) & (np.abs(r) < 1e10)
    r, p = r[mask], p[mask]
    if len(r) < 100:
        return None

    diff = r - p
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    bias = float(np.mean(diff))
    pearson = float(np.corrcoef(r.ravel(), p.ravel())[0, 1]) if np.std(r) > 0 and np.std(p) > 0 else 0.0
    pct5 = float(np.mean(np.abs(diff) < 0.05) * 100)
    pct1 = float(np.mean(np.abs(diff) < 0.01) * 100)
    return {
        "pearson_r": pearson, "rmse": rmse, "mae": mae, "bias": bias,
        "pct_within_0.05": pct5, "pct_within_0.01": pct1,
        "n_pixels": int(len(r)),
        "rust_mean": float(np.mean(r)), "py_mean": float(np.mean(p)),
    }


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def rust_binary():
    r = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io", "--example", "process_landsat"],
        cwd=str(REPO), capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        pytest.skip(f"Build failed: {r.stderr[:200]}")
    return REPO / "target" / "release" / "examples" / "process_landsat"


@pytest.fixture(scope="module")
def py_l8():
    return _run_python_full("L8")

@pytest.fixture(scope="module")
def py_l9():
    return _run_python_full("L9")

@pytest.fixture(scope="module")
def rust_l8(rust_binary):
    return _run_rust_full(rust_binary, "L8")

@pytest.fixture(scope="module")
def rust_l9(rust_binary):
    return _run_rust_full(rust_binary, "L9")


# ── Benchmark Tests ──────────────────────────────────────────────────────

class TestFullSceneBenchmark:
    """Time both pipelines on full scenes and report speedup."""

    def test_benchmark_l8(self, py_l8, rust_l8):
        py_t = py_l8["time"]
        rust_t = rust_l8["time"]
        if py_t is None:
            pytest.skip("No Python timing")
        speedup = py_t / rust_t if rust_t > 0 else 0
        print(f"\n  L8 Full Scene Benchmark:")
        print(f"    Python: {py_t:.1f}s")
        print(f"    Rust:   {rust_t:.1f}s")
        print(f"    Speedup: {speedup:.1f}x")
        if "timings" in rust_l8:
            t = rust_l8["timings"]
            print(f"    Rust breakdown: Load={t.get('load','?')}s AC={t.get('ac','?')}s Write={t.get('write','?')}s")

    def test_benchmark_l9(self, py_l9, rust_l9):
        py_t = py_l9["time"]
        rust_t = rust_l9["time"]
        if py_t is None:
            pytest.skip("No Python timing")
        speedup = py_t / rust_t if rust_t > 0 else 0
        print(f"\n  L9 Full Scene Benchmark:")
        print(f"    Python: {py_t:.1f}s")
        print(f"    Rust:   {rust_t:.1f}s")
        print(f"    Speedup: {speedup:.1f}x")
        if "timings" in rust_l9:
            t = rust_l9["timings"]
            print(f"    Rust breakdown: Load={t.get('load','?')}s AC={t.get('ac','?')}s Write={t.get('write','?')}s")

    def test_rust_faster_than_python_l8(self, py_l8, rust_l8):
        """Rust should be faster than Python on full scene."""
        if py_l8["time"] is None:
            pytest.skip("No Python timing")
        assert rust_l8["time"] < py_l8["time"], \
            f"Rust ({rust_l8['time']:.1f}s) not faster than Python ({py_l8['time']:.1f}s)"

    def test_rust_faster_than_python_l9(self, py_l9, rust_l9):
        if py_l9["time"] is None:
            pytest.skip("No Python timing")
        assert rust_l9["time"] < py_l9["time"], \
            f"Rust ({rust_l9['time']:.1f}s) not faster than Python ({py_l9['time']:.1f}s)"


# ── Full-Scene Accuracy ─────────────────────────────────────────────────

class TestFullSceneAccuracy:
    """Compare per-pixel rhos between Rust and Python on full scene."""

    @pytest.mark.parametrize("label", ["L8", "L9"])
    def test_all_bands_correlate(self, label, py_l8, py_l9, rust_l8, rust_l9, request):
        py = py_l8 if label == "L8" else py_l9
        rust = rust_l8 if label == "L8" else rust_l9
        cfg = SCENES[label]

        print(f"\n  {label} Full-Scene Accuracy:")
        print(f"  {'Band':<6} {'R':>10} {'RMSE':>10} {'Bias':>10} {'%<0.05':>8} {'%<0.01':>8}")
        print(f"  {'─'*6} {'─'*10} {'─'*10} {'─'*10} {'─'*8} {'─'*8}")

        all_r, all_rmse, all_pct5 = [], [], []
        for tif in rust["tifs"]:
            bname = tif.stem.split("_")[-1]
            wl = cfg["band_map"].get(bname)
            if not wl:
                continue
            py_data = _read_python_band(py["nc"], wl)
            if py_data is None:
                continue
            rust_data, rust_gt = _read_rust_band(tif)
            m = _align_and_compare(rust_data, rust_gt, py_data, py_data.shape)
            if m is None:
                continue
            print(f"  {bname:<6} {m['pearson_r']:>10.6f} {m['rmse']:>10.6f} "
                  f"{m['bias']:>10.6f} {m['pct_within_0.05']:>7.1f}% {m['pct_within_0.01']:>7.1f}%")
            all_r.append(m["pearson_r"])
            all_rmse.append(m["rmse"])
            all_pct5.append(m["pct_within_0.05"])

        assert len(all_r) >= 5, f"Only {len(all_r)} bands compared"
        mean_r = np.mean(all_r)
        mean_rmse = np.mean(all_rmse)
        mean_pct5 = np.mean(all_pct5)
        print(f"\n  Mean R={mean_r:.4f}  RMSE={mean_rmse:.4f}  %<0.05={mean_pct5:.1f}%")
        assert mean_r > 0.95, f"Mean R={mean_r:.4f} too low"


# ── Summary Report ───────────────────────────────────────────────────────

class TestBenchmarkSummary:
    def test_full_report(self, py_l8, py_l9, rust_l8, rust_l9):
        print("\n")
        print("=" * 78)
        print("  FULL-SCENE BENCHMARK: Rust vs Python ACOLITE (7931×7891 pixels)")
        print("  South Australia, Gulf St Vincent — Landsat 8 & 9")
        print("=" * 78)

        for label, py, rust in [("L8", py_l8, rust_l8), ("L9", py_l9, rust_l9)]:
            cfg = SCENES[label]
            py_t = py["time"]
            rust_t = rust["time"]
            speedup = py_t / rust_t if (py_t and rust_t and rust_t > 0) else None

            print(f"\n  {label}: {cfg['id']}")
            print(f"  {'─'*74}")
            print(f"    Python (full scene): {py_t:.1f}s" if py_t else "    Python: (cached)")
            print(f"    Rust   (full scene): {rust_t:.1f}s" if rust_t else "    Rust: (cached)")
            if speedup:
                print(f"    Speedup:             {speedup:.1f}x")
            if "timings" in rust:
                t = rust["timings"]
                print(f"    Rust breakdown:      Load={t.get('load','?')}s  "
                      f"AC={t.get('ac','?')}s  Write={t.get('write','?')}s")

        print(f"\n{'='*78}")
        print("  END BENCHMARK")
        print(f"{'='*78}")
