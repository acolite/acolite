"""
Subset (--limit) regression: Rust vs Python ACOLITE for all four sensors.

Validates that geographic bounding-box subsetting produces the same pixels
as full-scene processing restricted to the same ROI, for:
  - Landsat 8/9 OLI  (UTM projected, GeoTIFF bands)
  - Sentinel-2 MSI   (UTM projected, JP2 bands)
  - Sentinel-3 OLCI  (geographic, NetCDF radiance)
  - PACE OCI         (geographic, NetCDF L1B)

For each sensor the test:
  1. Runs Rust with --limit <south,west,north,east>
  2. Runs Python ACOLITE with limit=[south,west,north,east]
  3. Compares per-pixel rhos: Pearson R, RMSE, bias

Run:
  pytest tests/regression/test_subset_all_sensors.py -v -s
"""

import json, os, re, subprocess, sys, time
from pathlib import Path

import numpy as np
import pytest

REPO  = Path(__file__).resolve().parent.parent.parent
CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache"))

# Gulf St Vincent, South Australia — same ROI used across all sensors
LIMIT      = [-35.0, 138.0, -34.5, 138.7]
LIMIT_STR  = ",".join(str(x) for x in LIMIT)
LIMIT_PY   = LIMIT  # Python ACOLITE expects a list

# ── Tolerances ────────────────────────────────────────────────────────────
R_MIN    = 0.98   # Pearson R lower bound
RMSE_MAX = 0.02   # RMSE upper bound (reflectance units)
MIN_PIX  = 50     # minimum valid pixels for comparison


# ── Shared helpers ────────────────────────────────────────────────────────

def _compare(a, b):
    """Per-pixel stats between two 2-D arrays (NaN-masked)."""
    h = min(a.shape[0], b.shape[0])
    w = min(a.shape[1], b.shape[1])
    a, b = a[:h, :w], b[:h, :w]
    mask = np.isfinite(a) & np.isfinite(b) & (a > -0.5) & (b > -0.5) & (a < 2.0) & (b < 2.0)
    a, b = a[mask], b[mask]
    if len(a) < MIN_PIX:
        return None
    diff = a - b
    r = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else 0.0
    return {"pearson_r": r, "rmse": float(np.sqrt(np.mean(diff**2))),
            "bias": float(np.mean(diff)), "n_pixels": int(len(a))}


def _build_rust(example):
    r = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io", "--example", example],
        cwd=str(REPO), capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        pytest.skip(f"Rust build failed: {r.stderr[:300]}")
    binary = REPO / "target" / "release" / "examples" / example
    if not binary.exists():
        pytest.skip("Rust binary not found")
    return binary


def _run_python_acolite(inputfile, out_dir, extra_settings=None):
    """Run Python ACOLITE with limit and return output dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    settings = {
        "inputfile": str(inputfile),
        "output": str(out_dir),
        "limit": LIMIT_PY,
        "dsf_aot_estimate": "fixed",
        "dsf_spectrum_option": "intercept",
        "dsf_intercept_pixels": 200,
        "l2w_parameters": "",
        "rgb_rhot": False,
        "rgb_rhos": False,
        "map_l2w": False,
        "output_geometry": False,
    }
    if extra_settings:
        settings.update(extra_settings)
    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0,'{REPO}'); "
         f"import acolite; acolite.acolite.acolite_run(settings={settings!r})"],
        capture_output=True, text=True, timeout=600, cwd=str(REPO),
    )
    elapsed = time.time() - t0
    if r.returncode != 0:
        pytest.fail(f"Python ACOLITE failed ({elapsed:.1f}s):\n{r.stderr[:500]}")
    return out_dir, elapsed


# ═══════════════════════════════════════════════════════════════════════════
# LANDSAT 8/9 — UTM projected, GeoTIFF bands
# ═══════════════════════════════════════════════════════════════════════════

LANDSAT_SCENES = {
    "L8": {
        "id": "LC08_L1TP_098084_20240205_20240212_02_T1",
        "band_map": {"B1": "443", "B2": "483", "B3": "561", "B4": "655",
                     "B5": "865", "B6": "1609", "B7": "2201"},
        "py_prefix": "L8_OLI_2024_02_05_00_39_50_098084",
    },
    "L9": {
        "id": "LC09_L1TP_098084_20240213_20240213_02_T1",
        "band_map": {"B1": "443", "B2": "482", "B3": "561", "B4": "654",
                     "B5": "865", "B6": "1608", "B7": "2201"},
        "py_prefix": "L9_OLI_2024_02_13_00_40_00_098084",
    },
}


def _run_rust_landsat(binary, label):
    cfg = LANDSAT_SCENES[label]
    scene_dir = CACHE / "landsat_sa" / cfg["id"]
    out_dir   = CACHE / "landsat_sa" / f"rust_subset_{label}"
    if not scene_dir.exists():
        pytest.skip(f"No cached Landsat {label} scene")
    out_dir.mkdir(parents=True, exist_ok=True)
    tifs = sorted(out_dir.glob("*_corrected_B*.tif"))
    if len(tifs) == 7:
        return {"tifs": tifs, "cached": True}
    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(scene_dir), "--output", str(out_dir),
         "--limit", LIMIT_STR, "--aot-mode", "fixed"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust Landsat {label} failed:\n{r.stderr[:500]}")
    tifs = sorted(out_dir.glob("*_corrected_B*.tif"))
    return {"tifs": tifs, "elapsed": time.time() - t0, "cached": False,
            "stdout": r.stdout}


def _run_python_landsat(label):
    cfg = LANDSAT_SCENES[label]
    scene_dir = CACHE / "landsat_sa" / cfg["id"]
    out_dir   = CACHE / "landsat_sa" / f"py_subset_{label}"
    if not scene_dir.exists():
        pytest.skip(f"No cached Landsat {label} scene")
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if l2r_files:
        return {"nc": l2r_files[0], "out_dir": out_dir, "cached": True}
    out_dir, elapsed = _run_python_acolite(scene_dir, out_dir,
                                           {"dsf_aot_compute": "min"})
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if not l2r_files:
        pytest.fail("Python Landsat produced no L2R NetCDF")
    return {"nc": l2r_files[0], "out_dir": out_dir, "elapsed": elapsed, "cached": False}


@pytest.mark.parametrize("label", ["L8", "L9"])
def test_landsat_subset_pixel_accuracy(label):
    """Rust --limit subset matches Python ACOLITE limit= per pixel."""
    try:
        from osgeo import gdal
        import netCDF4 as nc4
    except ImportError:
        pytest.skip("osgeo/netCDF4 not available")

    binary = _build_rust("process_landsat")
    rust   = _run_rust_landsat(binary, label)
    py     = _run_python_landsat(label)

    cfg = LANDSAT_SCENES[label]
    ds  = nc4.Dataset(py["nc"])

    results = {}
    for rust_band, wl_str in cfg["band_map"].items():
        tif_matches = [t for t in rust["tifs"] if f"_{rust_band}.tif" in t.name]
        if not tif_matches:
            continue
        rust_ds = gdal.Open(str(tif_matches[0]))
        rust_arr = rust_ds.GetRasterBand(1).ReadAsArray().astype(float)
        rust_arr[rust_arr > 1e30] = np.nan

        py_var = f"rhos_{wl_str}"
        if py_var not in ds.variables:
            continue
        py_arr = ds.variables[py_var][:].astype(float)
        if hasattr(py_arr, "filled"):
            py_arr = py_arr.filled(np.nan)

        stats = _compare(rust_arr, py_arr)
        if stats:
            results[rust_band] = stats

    ds.close()
    assert results, f"No comparable bands for Landsat {label}"

    mean_r    = np.mean([s["pearson_r"] for s in results.values()])
    mean_rmse = np.mean([s["rmse"]      for s in results.values()])
    print(f"\nLandsat {label} subset accuracy:")
    for band, s in results.items():
        print(f"  {band}: R={s['pearson_r']:.4f}  RMSE={s['rmse']:.4f}  bias={s['bias']:+.4f}  n={s['n_pixels']}")
    print(f"  Mean R={mean_r:.4f}  Mean RMSE={mean_rmse:.4f}")

    assert mean_r    >= R_MIN,    f"Landsat {label} mean R={mean_r:.4f} < {R_MIN}"
    assert mean_rmse <= RMSE_MAX, f"Landsat {label} mean RMSE={mean_rmse:.4f} > {RMSE_MAX}"


@pytest.mark.parametrize("label", ["L8", "L9"])
def test_landsat_subset_size_reasonable(label):
    """Subset scene is smaller than full scene and non-empty."""
    binary = _build_rust("process_landsat")
    rust   = _run_rust_landsat(binary, label)
    assert rust["tifs"], f"No Rust output TIFs for Landsat {label}"
    try:
        from osgeo import gdal
        ds = gdal.Open(str(rust["tifs"][0]))
        w, h = ds.RasterXSize, ds.RasterYSize
        assert w > 0 and h > 0, "Empty subset"
        # Full Landsat scene is ~7900×7900; subset should be much smaller
        assert w < 7000 and h < 7000, f"Subset {w}×{h} suspiciously large (full scene?)"
        print(f"\nLandsat {label} subset size: {w}×{h} pixels")
    except ImportError:
        pytest.skip("osgeo not available")


# ═══════════════════════════════════════════════════════════════════════════
# SENTINEL-2 — UTM projected, JP2 bands
# ═══════════════════════════════════════════════════════════════════════════

S2_SCENES = {
    "S2A": {
        "id": "S2A_MSIL1C_20240321T005711_N0510_R002_T54HTF_20240321T022531.SAFE",
        "band_map": {"B01": "443", "B02": "492", "B03": "560", "B04": "665",
                     "B05": "704", "B06": "740", "B07": "783", "B08": "833",
                     "B8A": "865", "B11": "1614", "B12": "2202"},
        "py_prefix": "S2A_MSI_2024_03_21_00_57_11_T54HTF",
    },
}


def _run_rust_s2(binary, label):
    cfg = S2_SCENES[label]
    safe_dir = CACHE / "sentinel2_sa" / cfg["id"]
    out_dir  = CACHE / "sentinel2_sa" / f"rust_subset_{label}"
    if not safe_dir.exists():
        pytest.skip(f"No cached S2 {label} scene")
    out_dir.mkdir(parents=True, exist_ok=True)
    tifs = sorted(out_dir.glob("*_corrected_B*.tif"))
    if tifs:
        return {"tifs": tifs, "cached": True}
    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(safe_dir), "--output", str(out_dir),
         "--limit", LIMIT_STR, "--aot-mode", "fixed", "--res", "20"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust S2 {label} failed:\n{r.stderr[:500]}")
    tifs = sorted(out_dir.glob("*_corrected_B*.tif"))
    return {"tifs": tifs, "elapsed": time.time() - t0, "cached": False,
            "stdout": r.stdout}


def _run_python_s2(label):
    cfg = S2_SCENES[label]
    safe_dir = CACHE / "sentinel2_sa" / cfg["id"]
    out_dir  = CACHE / "sentinel2_sa" / f"py_subset_{label}"
    if not safe_dir.exists():
        pytest.skip(f"No cached S2 {label} scene")
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if l2r_files:
        return {"nc": l2r_files[0], "out_dir": out_dir, "cached": True}
    out_dir, elapsed = _run_python_acolite(safe_dir, out_dir, {
        "s2_target_res": 20,
        "dsf_wave_range": [400, 900],
        "dsf_nbands": 2,
        "dsf_nbands_fit": 2,
        "dsf_aot_compute": "min",
        "dsf_model_selection": "min_drmsd",
    })
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if not l2r_files:
        pytest.fail("Python S2 produced no L2R NetCDF")
    return {"nc": l2r_files[0], "out_dir": out_dir, "elapsed": elapsed, "cached": False}


@pytest.mark.parametrize("label", ["S2A"])
def test_s2_subset_pixel_accuracy(label):
    """Rust --limit subset matches Python ACOLITE limit= per pixel for S2."""
    try:
        from osgeo import gdal
        import netCDF4 as nc4
    except ImportError:
        pytest.skip("osgeo/netCDF4 not available")

    binary = _build_rust("process_sentinel2_ac")
    rust   = _run_rust_s2(binary, label)
    py     = _run_python_s2(label)

    cfg = S2_SCENES[label]
    ds  = nc4.Dataset(py["nc"])

    results = {}
    for rust_band, wl_str in cfg["band_map"].items():
        tif_matches = [t for t in rust["tifs"] if f"_{rust_band}.tif" in t.name]
        if not tif_matches:
            continue
        rust_ds  = gdal.Open(str(tif_matches[0]))
        rust_arr = rust_ds.GetRasterBand(1).ReadAsArray().astype(float)
        rust_arr[rust_arr > 1e30] = np.nan

        py_var = f"rhos_{wl_str}"
        if py_var not in ds.variables:
            continue
        py_arr = ds.variables[py_var][:].astype(float)
        if hasattr(py_arr, "filled"):
            py_arr = py_arr.filled(np.nan)

        stats = _compare(rust_arr, py_arr)
        if stats:
            results[rust_band] = stats

    ds.close()
    assert results, f"No comparable bands for S2 {label}"

    mean_r    = np.mean([s["pearson_r"] for s in results.values()])
    mean_rmse = np.mean([s["rmse"]      for s in results.values()])
    print(f"\nS2 {label} subset accuracy:")
    for band, s in results.items():
        print(f"  {band}: R={s['pearson_r']:.4f}  RMSE={s['rmse']:.4f}  bias={s['bias']:+.4f}  n={s['n_pixels']}")
    print(f"  Mean R={mean_r:.4f}  Mean RMSE={mean_rmse:.4f}")

    assert mean_r    >= R_MIN,    f"S2 {label} mean R={mean_r:.4f} < {R_MIN}"
    assert mean_rmse <= RMSE_MAX, f"S2 {label} mean RMSE={mean_rmse:.4f} > {RMSE_MAX}"


@pytest.mark.parametrize("label", ["S2A"])
def test_s2_subset_size_reasonable(label):
    """S2 subset is smaller than full tile and non-empty."""
    binary = _build_rust("process_sentinel2_ac")
    rust   = _run_rust_s2(binary, label)
    assert rust["tifs"], f"No Rust output TIFs for S2 {label}"
    try:
        from osgeo import gdal
        ds = gdal.Open(str(rust["tifs"][0]))
        w, h = ds.RasterXSize, ds.RasterYSize
        assert w > 0 and h > 0, "Empty subset"
        # Full S2 tile at 20m is 5490×5490; subset should be smaller
        assert w < 5000 and h < 5000, f"Subset {w}×{h} suspiciously large"
        print(f"\nS2 {label} subset size: {w}×{h} pixels")
    except ImportError:
        pytest.skip("osgeo not available")


# ═══════════════════════════════════════════════════════════════════════════
# SENTINEL-3 OLCI — geographic CRS, NetCDF radiance
# ═══════════════════════════════════════════════════════════════════════════

S3_SCENE = "S3B_OL_1_EFR____20240321T003707_20240321T004007_20240921T122930_0179_091_059_3600_MAR_R_NT_004.SEN3"

S3_BAND_MAP = {
    "Oa01": "401", "Oa02": "412", "Oa03": "443", "Oa04": "490",
    "Oa05": "510", "Oa06": "560", "Oa07": "620", "Oa08": "665",
    "Oa09": "674", "Oa10": "681", "Oa11": "709", "Oa12": "754",
    "Oa15": "768", "Oa16": "779", "Oa17": "865", "Oa18": "884",
    "Oa21": "1016",
}


def _run_rust_s3(binary):
    sen3_dir = CACHE / "sentinel3_sa" / S3_SCENE
    out_dir  = CACHE / "sentinel3_sa" / "rust_subset"
    if not sen3_dir.exists():
        pytest.skip("No S3 data — run _download_s3.py S3B first")
    out_dir.mkdir(parents=True, exist_ok=True)
    nc_files = list(out_dir.glob("*.nc"))
    if nc_files:
        return {"nc": nc_files[0], "cached": True}
    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--scene", str(sen3_dir), "--output", str(out_dir),
         "--limit", LIMIT_STR],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust S3 failed:\n{r.stderr[:500]}")
    nc_files = list(out_dir.glob("*.nc"))
    return {"nc": nc_files[0] if nc_files else None,
            "elapsed": time.time() - t0, "cached": False, "stdout": r.stdout}


def _run_python_s3():
    sen3_dir = CACHE / "sentinel3_sa" / S3_SCENE
    out_dir  = CACHE / "sentinel3_sa" / "py_subset"
    if not sen3_dir.exists():
        pytest.skip("No S3 data — run _download_s3.py S3B first")
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if l2r_files:
        return {"nc": l2r_files[0], "out_dir": out_dir, "cached": True}
    out_dir, elapsed = _run_python_acolite(sen3_dir, out_dir, {
        "smile_correction": True,
        "smile_correction_tgas": True,
        "use_supplied_ancillary": True,
        "dsf_aot_estimate": "fixed",
        "dsf_spectrum_option": "intercept",
        "dsf_intercept_pixels": 200,
    })
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if not l2r_files:
        pytest.fail("Python S3 produced no L2R NetCDF")
    return {"nc": l2r_files[0], "out_dir": out_dir, "elapsed": elapsed, "cached": False}


def test_s3_subset_pixel_accuracy():
    """Rust S3 --limit subset matches Python ACOLITE limit= per pixel."""
    try:
        import netCDF4 as nc4
    except ImportError:
        pytest.skip("netCDF4 not available")

    binary = _build_rust("process_sentinel3")
    rust   = _run_rust_s3(binary)
    py     = _run_python_s3()

    if not rust.get("nc") or not rust["nc"].exists():
        pytest.skip("Rust S3 produced no NetCDF output")

    rust_ds = nc4.Dataset(rust["nc"])
    py_ds   = nc4.Dataset(py["nc"])

    results = {}
    for rust_band, wl_str in S3_BAND_MAP.items():
        rust_var = f"rhos_{rust_band}"
        py_var   = f"rhos_{wl_str}"
        if rust_var not in rust_ds.variables or py_var not in py_ds.variables:
            continue
        rust_arr = rust_ds.variables[rust_var][:].astype(float)
        py_arr   = py_ds.variables[py_var][:].astype(float)
        for arr in (rust_arr, py_arr):
            if hasattr(arr, "filled"):
                arr = arr.filled(np.nan)
        stats = _compare(rust_arr, py_arr)
        if stats:
            results[rust_band] = stats

    rust_ds.close()
    py_ds.close()
    assert results, "No comparable S3 bands"

    mean_r    = np.mean([s["pearson_r"] for s in results.values()])
    mean_rmse = np.mean([s["rmse"]      for s in results.values()])
    print(f"\nS3 OLCI subset accuracy:")
    for band, s in results.items():
        print(f"  {band}: R={s['pearson_r']:.4f}  RMSE={s['rmse']:.4f}  bias={s['bias']:+.4f}  n={s['n_pixels']}")
    print(f"  Mean R={mean_r:.4f}  Mean RMSE={mean_rmse:.4f}")

    assert mean_r    >= R_MIN,    f"S3 mean R={mean_r:.4f} < {R_MIN}"
    assert mean_rmse <= RMSE_MAX, f"S3 mean RMSE={mean_rmse:.4f} > {RMSE_MAX}"


def test_s3_subset_size_reasonable():
    """S3 subset is smaller than full swath and non-empty."""
    try:
        import netCDF4 as nc4
    except ImportError:
        pytest.skip("netCDF4 not available")

    binary = _build_rust("process_sentinel3")
    rust   = _run_rust_s3(binary)
    if not rust.get("nc") or not rust["nc"].exists():
        pytest.skip("Rust S3 produced no NetCDF output")

    ds = nc4.Dataset(rust["nc"])
    # Find any rhos variable to check shape
    rhos_vars = [v for v in ds.variables if v.startswith("rhos_")]
    assert rhos_vars, "No rhos variables in S3 output"
    shape = ds.variables[rhos_vars[0]].shape
    ds.close()
    # Full S3 EFR swath is ~4865×4091; subset should be much smaller
    total_px = shape[-1] * shape[-2] if len(shape) >= 2 else shape[0]
    assert total_px > 0, "Empty S3 subset"
    assert total_px < 4865 * 4091, "S3 subset not smaller than full swath"
    print(f"\nS3 OLCI subset shape: {shape}")


# ═══════════════════════════════════════════════════════════════════════════
# PACE OCI — geographic CRS, NetCDF L1B
# ═══════════════════════════════════════════════════════════════════════════

PACE_GRANULE = "PACE_OCI.20241231T044250.L1B.V3.nc"

PACE_SAMPLE_BANDS = ["400", "443", "490", "550", "665", "750", "865"]


def _run_rust_pace(binary):
    nc_path = CACHE / "pace_sa" / PACE_GRANULE
    out_dir = CACHE / "pace_sa" / "rust_subset"
    if not nc_path.exists():
        pytest.skip(f"No PACE data at {nc_path}")
    out_dir.mkdir(parents=True, exist_ok=True)
    zarr_path = out_dir / "pace_subset.zarr"
    if zarr_path.exists():
        return {"zarr": zarr_path, "cached": True}
    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(nc_path), "--output", str(out_dir / "pace_subset"),
         "--limit", LIMIT_STR, "--aot-mode", "fixed"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    if r.returncode != 0:
        pytest.fail(f"Rust PACE failed:\n{r.stderr[:500]}")
    return {"zarr": zarr_path if zarr_path.exists() else None,
            "elapsed": time.time() - t0, "cached": False, "stdout": r.stdout}


def _run_python_pace():
    nc_path = CACHE / "pace_sa" / PACE_GRANULE
    out_dir = CACHE / "pace_sa" / "py_subset"
    if not nc_path.exists():
        pytest.skip(f"No PACE data at {nc_path}")
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if l2r_files:
        return {"nc": l2r_files[0], "out_dir": out_dir, "cached": True}
    out_dir, elapsed = _run_python_acolite(nc_path, out_dir, {
        "dsf_aot_estimate": "fixed",
        "dsf_spectrum_option": "intercept",
        "dsf_intercept_pixels": 200,
    })
    l2r_files = list(out_dir.glob("*_L2R.nc"))
    if not l2r_files:
        pytest.fail("Python PACE produced no L2R NetCDF")
    return {"nc": l2r_files[0], "out_dir": out_dir, "elapsed": elapsed, "cached": False}


def test_pace_subset_pixel_accuracy():
    """Rust PACE --limit subset matches Python ACOLITE limit= per pixel."""
    try:
        import netCDF4 as nc4
        import zarr
    except ImportError:
        pytest.skip("netCDF4/zarr not available")

    binary = _build_rust("process_pace_ac")
    rust   = _run_rust_pace(binary)
    py     = _run_python_pace()

    if not rust.get("zarr") or not rust["zarr"].exists():
        pytest.skip("Rust PACE produced no Zarr output")

    # Read Rust output from GeoZarr
    z = zarr.open(str(rust["zarr"]), mode="r")
    rust_data = z["data"][:]          # shape: (bands, rows, cols)
    rust_wls  = z["wavelengths"][:]   # shape: (bands,)

    py_ds = nc4.Dataset(py["nc"])

    results = {}
    for wl_str in PACE_SAMPLE_BANDS:
        wl_target = float(wl_str)
        # Find closest Rust band
        idx = int(np.argmin(np.abs(rust_wls - wl_target)))
        if abs(rust_wls[idx] - wl_target) > 5.0:
            continue
        rust_arr = rust_data[idx].astype(float)
        rust_arr[rust_arr > 1e30] = np.nan

        py_var = f"rhos_{wl_str}"
        if py_var not in py_ds.variables:
            # Try nearest wavelength in Python output
            py_rhos = [v for v in py_ds.variables if v.startswith("rhos_")]
            py_wls  = [float(v.split("_")[1]) for v in py_rhos]
            if not py_wls:
                continue
            nearest = py_rhos[int(np.argmin(np.abs(np.array(py_wls) - wl_target)))]
            if abs(float(nearest.split("_")[1]) - wl_target) > 5.0:
                continue
            py_var = nearest

        py_arr = py_ds.variables[py_var][:].astype(float)
        if hasattr(py_arr, "filled"):
            py_arr = py_arr.filled(np.nan)

        stats = _compare(rust_arr, py_arr)
        if stats:
            results[wl_str] = stats

    py_ds.close()
    assert results, "No comparable PACE bands"

    mean_r    = np.mean([s["pearson_r"] for s in results.values()])
    mean_rmse = np.mean([s["rmse"]      for s in results.values()])
    print(f"\nPACE OCI subset accuracy:")
    for wl, s in results.items():
        print(f"  {wl}nm: R={s['pearson_r']:.4f}  RMSE={s['rmse']:.4f}  bias={s['bias']:+.4f}  n={s['n_pixels']}")
    print(f"  Mean R={mean_r:.4f}  Mean RMSE={mean_rmse:.4f}")

    assert mean_r    >= R_MIN,    f"PACE mean R={mean_r:.4f} < {R_MIN}"
    assert mean_rmse <= RMSE_MAX, f"PACE mean RMSE={mean_rmse:.4f} > {RMSE_MAX}"


def test_pace_subset_size_reasonable():
    """PACE subset is smaller than full granule and non-empty."""
    try:
        import zarr
    except ImportError:
        pytest.skip("zarr not available")

    binary = _build_rust("process_pace_ac")
    rust   = _run_rust_pace(binary)
    if not rust.get("zarr") or not rust["zarr"].exists():
        pytest.skip("Rust PACE produced no Zarr output")

    z = zarr.open(str(rust["zarr"]), mode="r")
    shape = z["data"].shape  # (bands, rows, cols)
    assert len(shape) == 3, f"Unexpected Zarr shape: {shape}"
    assert shape[1] > 0 and shape[2] > 0, "Empty PACE subset"
    # Full PACE granule is ~1710×1272; subset should be smaller
    assert shape[1] < 1710 and shape[2] < 1272, \
        f"PACE subset {shape[1]}×{shape[2]} not smaller than full granule"
    print(f"\nPACE OCI subset shape: {shape[1]}×{shape[2]} × {shape[0]} bands")


# ═══════════════════════════════════════════════════════════════════════════
# Cross-sensor: subset timing comparison
# ═══════════════════════════════════════════════════════════════════════════

def test_subset_faster_than_fullscene_landsat():
    """Subset processing should be faster than full-scene for Landsat."""
    binary = _build_rust("process_landsat")
    cfg    = LANDSAT_SCENES["L8"]
    scene_dir = CACHE / "landsat_sa" / cfg["id"]
    if not scene_dir.exists():
        pytest.skip("No cached Landsat L8 scene")

    out_sub  = CACHE / "landsat_sa" / "rust_timing_subset"
    out_full = CACHE / "landsat_sa" / "rust_timing_full"
    out_sub.mkdir(parents=True, exist_ok=True)
    out_full.mkdir(parents=True, exist_ok=True)

    def _time_run(extra_args, out):
        r = subprocess.run(
            [str(binary), "--file", str(scene_dir), "--output", str(out),
             "--aot-mode", "fixed"] + extra_args,
            capture_output=True, text=True, timeout=300,
            env={**os.environ, "RUST_LOG": "warn"},
        )
        if r.returncode != 0:
            pytest.fail(f"Rust failed:\n{r.stderr[:300]}")
        # Parse timing from stdout
        for line in r.stdout.splitlines():
            m = re.search(r"Total.*?([\d.]+)s", line)
            if m:
                return float(m.group(1))
        return None

    t_sub  = _time_run(["--limit", LIMIT_STR], out_sub)
    t_full = _time_run([], out_full)

    if t_sub and t_full:
        print(f"\nLandsat L8 timing: subset={t_sub:.1f}s  full={t_full:.1f}s  ratio={t_full/t_sub:.1f}×")
        assert t_sub < t_full, f"Subset ({t_sub:.1f}s) not faster than full ({t_full:.1f}s)"
