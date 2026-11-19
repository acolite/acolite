"""
Full-scene benchmark: Rust vs Python ACOLITE for Sentinel-2 A/B.

Runs BOTH pipelines on complete S2 tiles (5490×5490 at 20m) and measures:
  - Wall-clock time for each pipeline
  - Per-pixel accuracy (Pearson R, RMSE, bias)
  - Speedup factor

South Australia, Gulf St Vincent — MGRS tile T54HTF.
"""

import os, re, subprocess, sys, time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent.parent
CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache")) / "sentinel2_sa"

# Rust band name → Python rhos wavelength name
BAND_MAP_S2A = {
    "B01": "443", "B02": "492", "B03": "560", "B04": "665",
    "B05": "704", "B06": "740", "B07": "783", "B08": "833",
    "B8A": "865", "B11": "1614", "B12": "2202",
}
BAND_MAP_S2B = {
    "B01": "442", "B02": "492", "B03": "559", "B04": "665",
    "B05": "704", "B06": "739", "B07": "780", "B08": "833",
    "B8A": "864", "B11": "1610", "B12": "2186",
}

SCENES = {
    "S2A": {
        "safe": "S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450.SAFE",
        "stem": "S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450",
        "band_map": BAND_MAP_S2A,
        "py_prefix": "S2A_MSI_2024_03_21_00_47_04_T54HTF",
    },
    "S2B": {
        "safe": "S2B_MSIL1C_20240228T004659_N0510_R102_T54HTF_20240228T020941.SAFE",
        "stem": "S2B_MSIL1C_20240228T004659_N0510_R102_T54HTF_20240228T020941",
        "band_map": BAND_MAP_S2B,
        "py_prefix": "S2B_MSI_2024_02_28_00_56_58_T54HTF",
    },
}


# ── Download S2 data from AWS ────────────────────────────────────────────

def _download_s2_scene(label):
    """Download S2 SAFE from AWS S3 if not cached."""
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    if safe_dir.exists() and (safe_dir / "GRANULE").exists():
        granules = list((safe_dir / "GRANULE").iterdir())
        if granules:
            jp2s = list((granules[0] / "IMG_DATA").glob("*.jp2"))
            if len(jp2s) >= 13:
                return safe_dir

    # Run download script
    download_script = REPO / "tests" / "regression" / "_download_s2.py"
    if not download_script.exists():
        pytest.skip(f"No S2 download script and no cached data for {label}")

    r = subprocess.run(
        [sys.executable, str(download_script), label],
        capture_output=True, text=True, timeout=600,
    )
    if r.returncode != 0:
        pytest.skip(f"S2 download failed: {r.stderr[:300]}")
    return safe_dir


# ── Python ACOLITE processing ────────────────────────────────────────────

def _run_python_full(label):
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    out_dir = CACHE / f"py_full_{label}"

    if not safe_dir.exists():
        _download_s2_scene(label)
    if not safe_dir.exists():
        pytest.skip(f"No cached {label} scene")

    l2r_nc = out_dir / f"{cfg['py_prefix']}_L2R.nc"
    if l2r_nc.exists():
        log = out_dir / "timing.txt"
        py_time = float(log.read_text().strip()) if log.exists() else None
        return {"nc": l2r_nc, "out_dir": out_dir, "time": py_time, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)

    settings = out_dir / "settings.txt"
    settings.write_text(
        f"inputfile={safe_dir}\n"
        f"output={out_dir}\n"
        "dsf_aot_estimate=fixed\n"
        "dsf_spectrum_option=intercept\n"
        "dsf_intercept_pixels=200\n"
        "dsf_interface_reflectance=False\n"
        "dsf_fixed_lut=ACOLITE-LUT-202110-MOD2\n"
        "l2w_parameters=\n"
        "rgb_rhot=False\n"
        "rgb_rhos=False\n"
        "map_l2w=False\n"
        "s2_target_res=20\n"
        "geometry_type=grids\n"
        "output_geometry=False\n"
        "resolved_geometry=False\n"
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

    (out_dir / "timing.txt").write_text(f"{elapsed:.2f}")

    l2r_nc = list(out_dir.glob("*_L2R.nc"))
    if not l2r_nc:
        pytest.fail(f"No L2R NetCDF produced for {label}:\n{r.stdout[-500:]}")

    return {"nc": l2r_nc[0], "out_dir": out_dir, "time": elapsed, "cached": False}


# ── Rust processing ──────────────────────────────────────────────────────

def _run_rust_full(binary, label):
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    out_dir = CACHE / f"rust_fixed_{label}"

    if not safe_dir.exists():
        pytest.skip(f"No cached {label} scene")

    stem = cfg["stem"]
    tifs = sorted(out_dir.glob(f"{stem}_corrected_B*.tif")) if out_dir.exists() else []
    log = out_dir / "timing.txt"
    if len(tifs) >= 11 and log.exists():
        rust_time = float(log.read_text().strip())
        return {"tifs": tifs, "out_dir": out_dir, "time": rust_time, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(safe_dir), "--output", str(out_dir), "--res", "20",
         "--aot-mode", "fixed", "--model", "auto"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0

    if r.returncode != 0:
        pytest.fail(f"Rust {label} failed ({elapsed:.1f}s):\n{r.stderr[:500]}")

    timings = {"total": elapsed}
    for line in r.stdout.splitlines():
        for key in ("Load", "AC", "Write"):
            m = re.search(rf"{key}:\s*([\d.]+)s", line)
            if m:
                timings[key.lower()] = float(m.group(1))

    (out_dir / "timing.txt").write_text(f"{elapsed:.2f}")

    tifs = sorted(out_dir.glob(f"{stem}_corrected_B*.tif"))
    return {"tifs": tifs, "out_dir": out_dir, "time": elapsed, "timings": timings, "cached": False}


# ── Pixel comparison ─────────────────────────────────────────────────────

def _read_python_band(nc_path, wave_name):
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
    from osgeo import gdal
    ds = gdal.Open(str(tif_path))
    data = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    ds = None
    return data


def _compare(rust_data, py_data):
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
    bias = float(np.mean(diff))
    pearson = float(np.corrcoef(r.ravel(), p.ravel())[0, 1]) if np.std(r) > 0 and np.std(p) > 0 else 0.0
    pct5 = float(np.mean(np.abs(diff) < 0.05) * 100)
    pct1 = float(np.mean(np.abs(diff) < 0.01) * 100)
    return {
        "pearson_r": pearson, "rmse": rmse, "bias": bias,
        "pct_within_0.05": pct5, "pct_within_0.01": pct1,
        "n_pixels": int(len(r)),
    }


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def rust_binary():
    r = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io", "--example", "process_sentinel2_ac"],
        cwd=str(REPO), capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        pytest.skip(f"Build failed: {r.stderr[:200]}")
    return REPO / "target" / "release" / "examples" / "process_sentinel2_ac"


@pytest.fixture(scope="module")
def py_s2a():
    return _run_python_full("S2A")

@pytest.fixture(scope="module")
def py_s2b():
    return _run_python_full("S2B")

@pytest.fixture(scope="module")
def rust_s2a(rust_binary):
    return _run_rust_full(rust_binary, "S2A")

@pytest.fixture(scope="module")
def rust_s2b(rust_binary):
    return _run_rust_full(rust_binary, "S2B")


# ── Benchmark Tests ──────────────────────────────────────────────────────

class TestS2Benchmark:
    """Time both pipelines on full S2 tiles and report speedup."""

    def test_benchmark_s2a(self, py_s2a, rust_s2a):
        py_t = py_s2a["time"]
        rust_t = rust_s2a["time"]
        if py_t is None:
            pytest.skip("No Python timing")
        speedup = py_t / rust_t if rust_t > 0 else 0
        print(f"\n  S2A Full Scene Benchmark:")
        print(f"    Python: {py_t:.1f}s")
        print(f"    Rust:   {rust_t:.1f}s")
        print(f"    Speedup: {speedup:.1f}x")
        if "timings" in rust_s2a:
            t = rust_s2a["timings"]
            print(f"    Rust breakdown: Load={t.get('load','?')}s AC={t.get('ac','?')}s Write={t.get('write','?')}s")

    def test_benchmark_s2b(self, py_s2b, rust_s2b):
        py_t = py_s2b["time"]
        rust_t = rust_s2b["time"]
        if py_t is None:
            pytest.skip("No Python timing")
        speedup = py_t / rust_t if rust_t > 0 else 0
        print(f"\n  S2B Full Scene Benchmark:")
        print(f"    Python: {py_t:.1f}s")
        print(f"    Rust:   {rust_t:.1f}s")
        print(f"    Speedup: {speedup:.1f}x")
        if "timings" in rust_s2b:
            t = rust_s2b["timings"]
            print(f"    Rust breakdown: Load={t.get('load','?')}s AC={t.get('ac','?')}s Write={t.get('write','?')}s")

    def test_rust_faster_s2a(self, py_s2a, rust_s2a):
        if py_s2a["time"] is None:
            pytest.skip("No Python timing")
        assert rust_s2a["time"] < py_s2a["time"]

    def test_rust_faster_s2b(self, py_s2b, rust_s2b):
        if py_s2b["time"] is None:
            pytest.skip("No Python timing")
        assert rust_s2b["time"] < py_s2b["time"]


# ── Accuracy ─────────────────────────────────────────────────────────────

class TestS2Accuracy:
    """Compare per-pixel rhos between Rust and Python."""

    @pytest.mark.parametrize("label", ["S2A", "S2B"])
    def test_all_bands_correlate(self, label, py_s2a, py_s2b, rust_s2a, rust_s2b):
        py = py_s2a if label == "S2A" else py_s2b
        rust = rust_s2a if label == "S2A" else rust_s2b
        cfg = SCENES[label]

        print(f"\n  {label} Full-Scene Accuracy (5490×5490 at 20m):")
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
            rust_data = _read_rust_band(tif)
            m = _compare(rust_data, py_data)
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
        assert mean_r > 0.90, f"Mean R={mean_r:.4f} too low"


# ── Summary Report ───────────────────────────────────────────────────────

class TestS2BenchmarkSummary:
    def test_full_report(self, py_s2a, py_s2b, rust_s2a, rust_s2b):
        print("\n")
        print("=" * 78)
        print("  FULL-SCENE BENCHMARK: Rust vs Python ACOLITE (5490×5490 at 20m)")
        print("  South Australia, Gulf St Vincent — Sentinel-2 A & B, MGRS T54HTF")
        print("=" * 78)

        for label, py, rust in [("S2A", py_s2a, rust_s2a), ("S2B", py_s2b, rust_s2b)]:
            cfg = SCENES[label]
            py_t = py["time"]
            rust_t = rust["time"]
            speedup = py_t / rust_t if (py_t and rust_t and rust_t > 0) else None

            print(f"\n  {label}: {cfg['safe']}")
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


# ── Tiled-vs-Tiled comparison ────────────────────────────────────────────

def _run_python_tiled(label):
    """Run Python ACOLITE with dsf_aot_estimate=tiled."""
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    out_dir = CACHE / f"py_tiled_{label}"

    if not safe_dir.exists():
        pytest.skip(f"No cached {label} scene")

    l2r_nc = list(out_dir.glob("*_L2R.nc")) if out_dir.exists() else []
    if l2r_nc:
        return {"nc": l2r_nc[0], "out_dir": out_dir, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)
    settings = out_dir / "settings.txt"
    settings.write_text(
        f"inputfile={safe_dir}\n"
        f"output={out_dir}\n"
        "dsf_aot_estimate=tiled\n"
        "dsf_spectrum_option=intercept\n"
        "dsf_intercept_pixels=200\n"
        "dsf_interface_reflectance=False\n"
        "l2w_parameters=\n"
        "rgb_rhot=False\n"
        "rgb_rhos=False\n"
        "map_l2w=False\n"
        "s2_target_res=20\n"
        "geometry_type=grids\n"
        "output_geometry=False\n"
        "resolved_geometry=False\n"
    )

    r = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0,'{REPO}'); "
         f"import acolite; acolite.acolite.acolite_run(settings='{settings}')"],
        capture_output=True, text=True, timeout=600, cwd=str(REPO),
    )
    if r.returncode != 0:
        pytest.fail(f"Python tiled {label} failed:\n{r.stderr[:500]}")

    l2r_nc = list(out_dir.glob("*_L2R.nc"))
    if not l2r_nc:
        pytest.fail(f"No L2R produced for tiled {label}")
    return {"nc": l2r_nc[0], "out_dir": out_dir, "cached": False}


def _run_rust_tiled(binary, label):
    """Run Rust with --aot-mode=tiled --model=auto."""
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    out_dir = CACHE / f"rust_tiled_{label}"

    if not safe_dir.exists():
        pytest.skip(f"No cached {label} scene")

    stem = cfg["stem"]
    tifs = sorted(out_dir.glob(f"{stem}_corrected_B*.tif")) if out_dir.exists() else []
    if len(tifs) >= 11:
        return {"tifs": tifs, "out_dir": out_dir, "cached": True}

    out_dir.mkdir(parents=True, exist_ok=True)
    r = subprocess.run(
        [str(binary), "--file", str(safe_dir), "--output", str(out_dir), "--res", "20",
         "--aot-mode", "tiled", "--model", "auto"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust tiled {label} failed:\n{r.stderr[:500]}")

    tifs = sorted(out_dir.glob(f"{stem}_corrected_B*.tif"))
    return {"tifs": tifs, "out_dir": out_dir, "cached": False}


class TestS2TiledMode:
    """Tiled-vs-tiled comparison: both Python and Rust use dsf_aot_estimate=tiled."""

    @pytest.mark.parametrize("label", ["S2A", "S2B"])
    def test_s2_tiled_mode(self, label, rust_binary):
        py = _run_python_tiled(label)
        rust = _run_rust_tiled(rust_binary, label)
        cfg = SCENES[label]

        print(f"\n  {label} Tiled-vs-Tiled (5490×5490 at 20m):")
        print(f"  {'Band':<6} {'R':>10} {'RMSE':>10} {'Bias':>10} {'%<0.05':>8}")
        print(f"  {'─'*6} {'─'*10} {'─'*10} {'─'*10} {'─'*8}")

        all_r = []
        for tif in rust["tifs"]:
            bname = tif.stem.split("_")[-1]
            wl = cfg["band_map"].get(bname)
            if not wl:
                continue
            py_data = _read_python_band(py["nc"], wl)
            if py_data is None:
                continue
            rust_data = _read_rust_band(tif)
            m = _compare(rust_data, py_data)
            if m is None:
                continue
            print(f"  {bname:<6} {m['pearson_r']:>10.6f} {m['rmse']:>10.6f} "
                  f"{m['bias']:>10.6f} {m['pct_within_0.05']:>7.1f}%")
            all_r.append(m["pearson_r"])

        assert len(all_r) >= 5, f"Only {len(all_r)} bands compared"
        assert np.mean(all_r) > 0.90, f"Mean R={np.mean(all_r):.4f} too low"
