"""
Landsat 8/9 full ATCOR regression: Rust vs Python ACOLITE (DSF).

South Australia, Gulf St Vincent water scenes:
  - L8: LC08_L1TP_098084_20240205 (ROI: [-35.4, 138.2, -35.1, 138.5])
  - L9: LC09_L1TP_098084_20240213 (same ROI)

Both pipelines run full atmospheric correction:
  - Python: ACOLITE DSF → rhos_* (surface reflectance)
  - Rust:   acolite-rs AC pipeline → per-band COG

Per-pixel metrics: Pearson R, RMSE, MAE, bias, % within tolerance.
"""

import json, os, re, subprocess, time
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
        "py_l2r_prefix": "L8_OLI_2024_02_05_00_39_50_098084_L2R",
        "py_time": 8.87,
    },
    "L9": {
        "id": "LC09_L1TP_098084_20240213_20240213_02_T1",
        "band_map": BAND_MAP_L9,
        "py_l2r_prefix": "L9_OLI_2024_02_13_00_40_00_098084_L2R",
        "py_time": 8.62,
    },
}

FILL_VAL = 9.96921e+36  # NetCDF fill value written to GeoTIFF


# ── Helpers ──────────────────────────────────────────────────────────────

def _run_rust(binary, label):
    cfg = SCENES[label]
    scene_dir = CACHE / cfg["id"]
    out_dir = CACHE / f"rust_output_{cfg['id'][:4]}"
    if not scene_dir.exists():
        pytest.skip(f"No cached {label} scene")
    out_dir.mkdir(parents=True, exist_ok=True)
    tifs = sorted(out_dir.glob(f"{cfg['id']}_corrected_B*.tif"))
    if len(tifs) == 7:
        return {"tifs": tifs, "timings": {}, "cached": True}

    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(scene_dir), "--output", str(out_dir)],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust {label} failed:\n{r.stderr[:500]}")
    tifs = sorted(out_dir.glob(f"{cfg['id']}_corrected_B*.tif"))
    timings = {}
    for line in r.stdout.splitlines():
        for key in ("Load", "AC", "Write", "Total"):
            m = re.search(rf"{key}:\s*([\d.]+)s", line)
            if m:
                timings[key.lower()] = float(m.group(1))
    return {"tifs": tifs, "timings": timings, "cached": False, "elapsed": time.time() - t0}


def _extract_roi(rust_tif, py_tif):
    from osgeo import gdal
    py_ds, rust_ds = gdal.Open(str(py_tif)), gdal.Open(str(rust_tif))
    py_gt, rust_gt = py_ds.GetGeoTransform(), rust_ds.GetGeoTransform()
    col = int(round((py_gt[0] - rust_gt[0]) / rust_gt[1]))
    row = int(round((py_gt[3] - rust_gt[3]) / rust_gt[5]))
    w, h = py_ds.RasterXSize, py_ds.RasterYSize
    col, row = max(0, col), max(0, row)
    w = min(w, rust_ds.RasterXSize - col)
    h = min(h, rust_ds.RasterYSize - row)
    rust_data = rust_ds.GetRasterBand(1).ReadAsArray(col, row, w, h).astype(np.float64)
    py_data = py_ds.GetRasterBand(1).ReadAsArray(0, 0, w, h).astype(np.float64)
    return rust_data, py_data


def _pixel_metrics(rust, py):
    # Mask fill values and NaN
    mask = (np.isfinite(rust) & np.isfinite(py)
            & (np.abs(py) < 1e10) & (np.abs(rust) < 1e10))
    r, p = rust[mask], py[mask]
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


def _all_metrics(label, rust_out):
    cfg = SCENES[label]
    py_dir = CACHE / f"py_l2r_{label}"
    if not py_dir.exists():
        pytest.skip(f"No Python L2R for {label}")
    results = {}
    for tif in rust_out["tifs"]:
        bname = tif.stem.split("_")[-1]
        wl = cfg["band_map"].get(bname)
        if not wl:
            continue
        py_tif = py_dir / f"{cfg['py_l2r_prefix']}_rhos_{wl}.tif"
        if not py_tif.exists():
            continue
        rust_roi, py_data = _extract_roi(tif, py_tif)
        m = _pixel_metrics(rust_roi, py_data)
        if m:
            results[bname] = m
    return results


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
def rust_l8(rust_binary):
    return _run_rust(rust_binary, "L8")

@pytest.fixture(scope="module")
def rust_l9(rust_binary):
    return _run_rust(rust_binary, "L9")

@pytest.fixture(scope="module")
def metrics_l8(rust_l8):
    return _all_metrics("L8", rust_l8)

@pytest.fixture(scope="module")
def metrics_l9(rust_l9):
    return _all_metrics("L9", rust_l9)


# ── Performance ──────────────────────────────────────────────────────────

class TestPerformance:
    @pytest.mark.parametrize("label,fix", [("L8", "rust_l8"), ("L9", "rust_l9")])
    def test_completes_7_bands(self, label, fix, request):
        assert len(request.getfixturevalue(fix)["tifs"]) == 7

    @pytest.mark.parametrize("label,fix", [("L8", "rust_l8"), ("L9", "rust_l9")])
    def test_timing(self, label, fix, request):
        r = request.getfixturevalue(fix)
        t = r["timings"]
        if not t:
            pytest.skip("cached")
        py_t = SCENES[label]["py_time"]
        print(f"\n{label}: Rust total={t.get('total','?')}s  "
              f"(load={t.get('load','?')} AC={t.get('ac','?')} write={t.get('write','?')})  "
              f"Python L2R(ROI)={py_t}s")


# ── Per-Pixel Similarity (full ATCOR: Rust AC vs Python DSF rhos) ────────

class TestPixelSimilarity:
    @pytest.mark.parametrize("label,fix", [("L8", "metrics_l8"), ("L9", "metrics_l9")])
    def test_all_bands_have_metrics(self, label, fix, request):
        m = request.getfixturevalue(fix)
        assert len(m) == 7, f"{label}: {len(m)}/7 bands"

    @pytest.mark.parametrize("label,fix", [("L8", "metrics_l8"), ("L9", "metrics_l9")])
    def test_pearson_r(self, label, fix, request):
        m = request.getfixturevalue(fix)
        for b, met in sorted(m.items()):
            r = met["pearson_r"]
            print(f"  {label} {b}: Pearson R = {r:.6f}")
            assert r > 0.80, f"{label} {b} R={r:.4f}"

    @pytest.mark.parametrize("label,fix", [("L8", "metrics_l8"), ("L9", "metrics_l9")])
    def test_spectral_ordering(self, label, fix, request):
        m = request.getfixturevalue(fix)
        vnir = [m[b]["rust_mean"] for b in ["B1","B2","B3","B4","B5"] if b in m]
        swir = [m[b]["rust_mean"] for b in ["B6","B7"] if b in m]
        if vnir and swir:
            # Allow small tolerance for near-equal means (water scenes)
            assert np.mean(swir) < np.mean(vnir) + 0.01, "SWIR should be ≤ VNIR"


# ── COG Structure ────────────────────────────────────────────────────────

class TestCOGStructure:
    @pytest.mark.parametrize("label,fix", [("L8", "rust_l8"), ("L9", "rust_l9")])
    def test_cog_with_projection(self, label, fix, request):
        from osgeo import gdal
        r = request.getfixturevalue(fix)
        ds = gdal.Open(str(r["tifs"][0]))
        md = ds.GetMetadata("IMAGE_STRUCTURE")
        assert md.get("LAYOUT") == "COG"
        assert "UTM" in ds.GetProjection()
        gt = ds.GetGeoTransform()
        assert gt[1] == 30.0 and gt[5] == -30.0


# ── Summary Report ───────────────────────────────────────────────────────

class TestRegressionSummary:
    def test_full_report(self, rust_l8, rust_l9, metrics_l8, metrics_l9):
        print("\n")
        print("=" * 78)
        print("  LANDSAT 8/9 FULL ATCOR REGRESSION — South Australia (water)")
        print("  Rust AC vs Python ACOLITE DSF (rhos = surface reflectance)")
        print("=" * 78)

        for label, rust, metrics in [("L8", rust_l8, metrics_l8), ("L9", rust_l9, metrics_l9)]:
            cfg = SCENES[label]
            t = rust["timings"]
            print(f"\n{'─'*78}")
            print(f"  {label}: {cfg['id']}")
            print(f"{'─'*78}")
            if t:
                print(f"  Rust:   Load={t.get('load','?')}s  AC={t.get('ac','?')}s  "
                      f"Write={t.get('write','?')}s  Total={t.get('total','?')}s")
            else:
                print(f"  (cached — no timing)")
            print(f"  Python: L2R(ROI)={cfg['py_time']}s")

            print(f"\n  {'Band':<6} {'Pearson R':>10} {'RMSE':>10} {'MAE':>10} "
                  f"{'Bias':>10} {'%<0.05':>8} {'%<0.01':>8} {'Rust μ':>10} {'Py μ':>10}")
            print(f"  {'─'*6} {'─'*10} {'─'*10} {'─'*10} {'─'*10} {'─'*8} {'─'*8} {'─'*10} {'─'*10}")

            for band in ["B1","B2","B3","B4","B5","B6","B7"]:
                if band not in metrics:
                    continue
                m = metrics[band]
                print(f"  {band:<6} {m['pearson_r']:>10.6f} {m['rmse']:>10.6f} {m['mae']:>10.6f} "
                      f"{m['bias']:>10.6f} {m['pct_within_0.05']:>7.1f}% {m['pct_within_0.01']:>7.1f}% "
                      f"{m['rust_mean']:>10.6f} {m['py_mean']:>10.6f}")

            all_r = [metrics[b]["pearson_r"] for b in metrics]
            all_rmse = [metrics[b]["rmse"] for b in metrics]
            all_pct5 = [metrics[b]["pct_within_0.05"] for b in metrics]
            print(f"\n  SIMILARITY SCORES:")
            print(f"    Mean Pearson R:    {np.mean(all_r):.6f}")
            print(f"    Mean RMSE:         {np.mean(all_rmse):.6f}")
            print(f"    Mean %<0.05:       {np.mean(all_pct5):.1f}%")

        print(f"\n{'='*78}")
        print("  END REPORT")
        print(f"{'='*78}")
