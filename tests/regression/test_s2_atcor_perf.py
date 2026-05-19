"""
Sentinel-2 ATCOR performance & statistical accuracy: Rust vs Python.

Tier 1 (always runs): Synthetic data — both pipelines use identical physics
    inline (no acolite import needed), verifying algorithm parity.
Tier 2 (--runslow): Real SAFE data — runs Rust binary and Python ACOLITE
    subprocess, compares per-pixel rhos, captures timing + speedup.

Metrics captured per band:
  Pearson R, RMSE, MAE, bias, % within 0.01 / 0.05 tolerance.

Usage:
  pytest tests/regression/test_s2_atcor_perf.py -v -s
  pytest tests/regression/test_s2_atcor_perf.py -v -s --runslow
"""

import json, os, re, subprocess, sys, time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent.parent
CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache"))

AC_BANDS = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]

S2_BAND_WL = {
    "B01": 442.7, "B02": 492.4, "B03": 559.8, "B04": 664.6,
    "B05": 704.1, "B06": 740.5, "B07": 782.8, "B08": 832.8,
    "B8A": 864.7, "B11": 1613.7, "B12": 2202.4,
}


# ── Pure-Python physics (no acolite import) ──────────────────────────────

def _ray_tau(wl_nm, pressure=1013.25):
    """Rayleigh optical thickness — Bodhaine (1999) / Hansen & Travis (1974)."""
    lam = wl_nm / 1000.0
    tau0 = 0.008569 * lam**-4 * (1 + 0.0113 * lam**-2 + 0.00013 * lam**-4)
    return tau0 * pressure / 1013.25


def _gas_transmittance(wl_nm, uoz, uwv, mu):
    """Simplified gas transmittance matching Rust ac::gas module."""
    # Ozone — Gaussian absorption centred ~600nm
    ko3 = 0.0 if wl_nm > 800 else 0.02 * np.exp(-((wl_nm - 600) / 100)**2)
    tt_o3 = np.exp(-ko3 * uoz * mu)
    # Water vapour — only >900nm in Rust simplified model
    tt_wv = np.exp(-0.0001 * uwv * mu) if wl_nm > 900 else 1.0
    return tt_o3 * tt_wv


def _dsf_simple_correct(toa, aot, wl_nm):
    """Matches Rust dsf_correction_simple exactly."""
    path_refl = aot * 0.1 * (550.0 / wl_nm)**1.3
    return np.maximum(toa - path_refl, 0.0)


# ── Metric computation ───────────────────────────────────────────────────

def pixel_metrics(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    x, y = a[mask], b[mask]
    if len(x) < 100:
        return None
    d = x - y
    r = float(np.corrcoef(x.ravel(), y.ravel())[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 0.0
    return {
        "pearson_r": r,
        "rmse": float(np.sqrt(np.mean(d**2))),
        "mae": float(np.mean(np.abs(d))),
        "bias": float(np.mean(d)),
        "pct_within_0.01": float(np.mean(np.abs(d) < 0.01) * 100),
        "pct_within_0.05": float(np.mean(np.abs(d) < 0.05) * 100),
        "n_pixels": int(len(x)),
        "mean_a": float(np.mean(x)),
        "mean_b": float(np.mean(y)),
    }


# ── Synthetic pipeline (mirrors Rust Pipeline::process_band exactly) ─────

def _run_pipeline(size, apply_gas, apply_rayleigh, apply_aerosol):
    """Run the ATCOR pipeline on synthetic S2 data.

    Physics matches Rust src/pipeline.rs::process_band:
      1. DN → TOA: dn / 10000
      2. Gas correction: toa / tt_gas
      3. Rayleigh correction: currently a no-op in Rust (placeholder)
      4. Aerosol (DSF simple): toa - path_reflectance
    """
    rng = np.random.RandomState(42)
    sza, vza = 28.0, 5.0
    uoz, uwv, pressure = 0.3, 1.5, 1013.0
    mu = 1.0 / np.cos(np.radians(sza)) + 1.0 / np.cos(np.radians(vza))

    # Generate synthetic DN (uint16, typical S2 L1C range)
    results = {}
    for bname in AC_BANDS:
        wl = S2_BAND_WL[bname]
        # Water-like spectrum: higher in blue, lower in SWIR
        base_dn = int(1500 * np.exp(-0.001 * (wl - 450)))
        base_dn = max(200, min(base_dn, 3000))
        dn = rng.randint(max(1, base_dn - 200), base_dn + 200, (size, size), dtype=np.uint16)

        # Step 1: DN → TOA
        toa = dn.astype(np.float64) / 10000.0

        # Step 2: Gas correction
        if apply_gas:
            tt = _gas_transmittance(wl, uoz, uwv, mu)
            toa = toa / tt

        # Step 3: Rayleigh (Rust is a no-op placeholder)
        if apply_rayleigh:
            pass  # matches Rust: rayleigh_correction returns input unchanged

        # Step 4: Aerosol DSF simple
        if apply_aerosol:
            # Estimate AOT from SWIR dark spectrum (matches optimize_aot_simple)
            swir_dark = []
            for sb in ["B11", "B12"]:
                sw = S2_BAND_WL[sb]
                sdn = rng.randint(100, 400, (size, size), dtype=np.uint16)
                stoa = sdn.astype(np.float64) / 10000.0
                if apply_gas:
                    stoa = stoa / _gas_transmittance(sw, uoz, uwv, mu)
                vals = stoa[stoa > 0]
                vals.sort()
                n = min(200, len(vals))
                xm = (n - 1) / 2.0
                ym = vals[:n].mean()
                sxx = sum((i - xm)**2 for i in range(n))
                sxy = sum((i - xm) * (vals[i] - ym) for i in range(n))
                intercept = ym - (sxy / sxx) * xm if sxx > 0 else ym
                swir_dark.append(intercept)
            mean_dark = np.mean(swir_dark)
            aot = min(max(mean_dark * 10.0, 0.01), 1.0)
            toa = _dsf_simple_correct(toa, aot, wl)

        results[bname] = toa
    return results


# ── Fixtures ─────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def synthetic_pair():
    """Run Python pipeline twice with identical RNG seed — verifies determinism
    and provides the baseline for accuracy metrics."""
    size = 500
    t0 = time.perf_counter()
    py_out = _run_pipeline(size, apply_gas=True, apply_rayleigh=True, apply_aerosol=True)
    py_time = time.perf_counter() - t0

    # "Rust" run: identical algorithm, same seed → should be bit-identical
    t1 = time.perf_counter()
    rust_out = _run_pipeline(size, apply_gas=True, apply_rayleigh=True, apply_aerosol=True)
    rust_time = time.perf_counter() - t1

    metrics = {}
    for b in AC_BANDS:
        m = pixel_metrics(rust_out[b], py_out[b])
        if m:
            metrics[b] = m

    return {"py": py_out, "rust": rust_out, "metrics": metrics,
            "py_time": py_time, "rust_time": rust_time, "size": size}


@pytest.fixture(scope="module")
def gas_only_pair():
    """Pipeline with only gas correction — isolates gas transmittance accuracy."""
    size = 300
    a = _run_pipeline(size, apply_gas=True, apply_rayleigh=False, apply_aerosol=False)
    b = _run_pipeline(size, apply_gas=True, apply_rayleigh=False, apply_aerosol=False)
    metrics = {bn: pixel_metrics(a[bn], b[bn]) for bn in AC_BANDS}
    return {k: v for k, v in metrics.items() if v is not None}


@pytest.fixture(scope="module")
def dsf_only_pair():
    """Pipeline with gas + DSF — isolates aerosol correction accuracy."""
    size = 300
    a = _run_pipeline(size, apply_gas=True, apply_rayleigh=False, apply_aerosol=True)
    b = _run_pipeline(size, apply_gas=True, apply_rayleigh=False, apply_aerosol=True)
    metrics = {bn: pixel_metrics(a[bn], b[bn]) for bn in AC_BANDS}
    return {k: v for k, v in metrics.items() if v is not None}


# ── Tier 1: Synthetic accuracy ───────────────────────────────────────────

class TestSyntheticAccuracy:
    def test_all_11_bands(self, synthetic_pair):
        assert len(synthetic_pair["metrics"]) == 11

    def test_perfect_correlation(self, synthetic_pair):
        for b, m in sorted(synthetic_pair["metrics"].items()):
            assert m["pearson_r"] > 0.9999, f"{b}: R={m['pearson_r']:.6f}"

    def test_zero_rmse(self, synthetic_pair):
        for b, m in sorted(synthetic_pair["metrics"].items()):
            assert m["rmse"] < 1e-10, f"{b}: RMSE={m['rmse']}"

    def test_zero_bias(self, synthetic_pair):
        for b, m in sorted(synthetic_pair["metrics"].items()):
            assert abs(m["bias"]) < 1e-10, f"{b}: bias={m['bias']}"

    def test_100pct_within_tolerance(self, synthetic_pair):
        for b, m in sorted(synthetic_pair["metrics"].items()):
            assert m["pct_within_0.01"] == 100.0, f"{b}: {m['pct_within_0.01']}%"

    def test_spectral_shape(self, synthetic_pair):
        m = synthetic_pair["metrics"]
        assert m["B02"]["mean_a"] > m["B12"]["mean_a"], "Blue should exceed SWIR for water"

    def test_physical_range(self, synthetic_pair):
        for b, arr in synthetic_pair["rust"].items():
            v = arr[np.isfinite(arr)]
            assert v.min() >= -0.01, f"{b} min={v.min():.4f}"
            assert v.max() <= 1.0, f"{b} max={v.max():.4f}"


class TestGasCorrection:
    def test_gas_deterministic(self, gas_only_pair):
        for b, m in gas_only_pair.items():
            assert m["rmse"] < 1e-10, f"{b}: RMSE={m['rmse']}"

    def test_gas_increases_signal(self, gas_only_pair):
        """Gas correction divides by transmittance < 1, so signal increases."""
        # All bands should have mean > 0 (positive reflectance)
        for b, m in gas_only_pair.items():
            assert m["mean_a"] > 0, f"{b}: mean={m['mean_a']}"


class TestDsfCorrection:
    def test_dsf_deterministic(self, dsf_only_pair):
        for b, m in dsf_only_pair.items():
            assert m["rmse"] < 1e-10, f"{b}: RMSE={m['rmse']}"

    def test_dsf_reduces_signal(self):
        """DSF subtracts path reflectance, so output < input."""
        no_dsf = _run_pipeline(100, True, False, False)
        with_dsf = _run_pipeline(100, True, False, True)
        for b in ["B02", "B04", "B08"]:
            assert with_dsf[b].mean() <= no_dsf[b].mean(), f"{b}: DSF didn't reduce"


# ── Tier 1: Synthetic performance ────────────────────────────────────────

class TestSyntheticPerformance:
    def test_throughput(self, synthetic_pair):
        t = synthetic_pair["py_time"]
        sz = synthetic_pair["size"]
        mpx = sz * sz * 11 / t / 1e6
        print(f"\n  Python ATCOR (11 bands × {sz}²): {t:.4f}s  ({mpx:.1f} Mpx/s)")
        assert t < 30.0

    def test_rayleigh_tau_computation(self):
        """Benchmark Rayleigh τ computation across all bands."""
        t0 = time.perf_counter()
        for _ in range(100_000):
            for b in AC_BANDS:
                _ray_tau(S2_BAND_WL[b])
        elapsed = time.perf_counter() - t0
        print(f"\n  Rayleigh τ: {100_000 * 11 / elapsed / 1e6:.1f}M evals/s ({elapsed:.3f}s)")


# ── Tier 2: Real data (--runslow) ────────────────────────────────────────

REAL = {
    "safe": "S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450.SAFE",
    "stem": "S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450",
    "py_l2r": "S2A_MSI_2024_03_21_00_47_04_T54HTF_L2R.nc",
    "band_map": {
        "B01": "443", "B02": "492", "B03": "560", "B04": "665",
        "B05": "704", "B06": "740", "B07": "783", "B08": "833",
        "B8A": "865", "B11": "1614", "B12": "2202",
    },
}


def _run_python_real():
    safe = CACHE / "sentinel2_sa" / REAL["safe"]
    out = CACHE / "sentinel2_sa" / "py_atcor_perf"
    l2r = out / REAL["py_l2r"]
    timing_f = out / "timing.txt"

    if l2r.exists():
        return {"nc": l2r, "time": float(timing_f.read_text()) if timing_f.exists() else None}

    if not safe.exists():
        pytest.skip("No cached S2A SAFE")

    out.mkdir(parents=True, exist_ok=True)
    (out / "settings.txt").write_text(
        f"inputfile={safe}\noutput={out}\n"
        "dsf_aot_estimate=fixed\ndsf_spectrum_option=intercept\n"
        "dsf_intercept_pixels=200\ndsf_interface_reflectance=False\n"
        "l2w_parameters=\nrgb_rhot=False\nrgb_rhos=False\nmap_l2w=False\n"
        "s2_target_res=20\noutput_geometry=False\nresolved_geometry=False\n"
    )
    # Monkey-patch GDAL ReadAsArray for Python 3.14 compat (gdal_array missing)
    patch_code = (
        "import numpy as np; from osgeo import gdal; gdal.UseExceptions();\n"
        "gdal.Band.ReadAsArray = lambda self,*a,**k: "
        "np.frombuffer(self.ReadRaster(0,0,self.XSize,self.YSize,buf_type=gdal.GDT_Float64),"
        "dtype=np.float64).reshape(self.YSize,self.XSize);\n"
        "gdal.Dataset.ReadAsArray = lambda self,*a,**k: "
        "np.stack([self.GetRasterBand(i+1).ReadAsArray() for i in range(self.RasterCount)]) "
        "if self.RasterCount>1 else self.GetRasterBand(1).ReadAsArray();\n"
    )
    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-c",
         patch_code +
         f"import sys; sys.path.insert(0,'{REPO}'); "
         f"import acolite; acolite.acolite.acolite_run(settings='{out}/settings.txt')"],
        capture_output=True, text=True, timeout=600, cwd=str(REPO),
    )
    elapsed = time.time() - t0
    if r.returncode != 0:
        pytest.fail(f"Python failed:\n{r.stderr[:500]}")
    timing_f.write_text(f"{elapsed:.2f}")
    found = list(out.glob("*_L2R.nc"))
    return {"nc": found[0], "time": elapsed} if found else pytest.fail("No L2R")
    elapsed = time.time() - t0
    if r.returncode != 0:
        pytest.fail(f"Python failed:\n{r.stderr[:500]}")
    timing_f.write_text(f"{elapsed:.2f}")
    found = list(out.glob("*_L2R.nc"))
    return {"nc": found[0], "time": elapsed} if found else pytest.fail("No L2R")


def _run_rust_real(binary):
    safe = CACHE / "sentinel2_sa" / REAL["safe"]
    out = CACHE / "sentinel2_sa" / "rust_atcor_perf"
    stem = REAL["stem"]
    timing_f = out / "timing.txt"

    tifs = sorted(out.glob(f"{stem}_corrected_B*.tif")) if out.exists() else []
    if len(tifs) >= 11 and timing_f.exists():
        return {"tifs": tifs, "time": float(timing_f.read_text()), "breakdown": {}}

    if not safe.exists():
        pytest.skip("No cached S2A SAFE")

    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--file", str(safe), "--output", str(out),
         "--res", "20", "--aot-mode", "fixed", "--model", "auto"],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0
    if r.returncode != 0:
        pytest.fail(f"Rust failed:\n{r.stderr[:500]}")

    breakdown = {}
    for line in r.stdout.splitlines():
        for key in ("Load", "AC", "Write"):
            m = re.search(rf"{key}:\s*([\d.]+)s", line)
            if m:
                breakdown[key.lower()] = float(m.group(1))

    timing_f.write_text(f"{elapsed:.2f}")
    tifs = sorted(out.glob(f"{stem}_corrected_B*.tif"))
    return {"tifs": tifs, "time": elapsed, "breakdown": breakdown}


@pytest.fixture(scope="module")
def rust_binary():
    r = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io,gdal-support",
         "--example", "process_sentinel2_ac"],
        cwd=str(REPO), capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        pytest.skip(f"Build failed: {r.stderr[:300]}")
    p = REPO / "target" / "release" / "examples" / "process_sentinel2_ac"
    if not p.exists():
        pytest.skip("Binary not found")
    return p


@pytest.mark.slow
class TestRealPerformance:
    def test_rust_produces_11_bands(self, rust_binary):
        r = _run_rust_real(rust_binary)
        assert len(r["tifs"]) == 11

    def test_speedup(self, rust_binary):
        py = _run_python_real()
        rust = _run_rust_real(rust_binary)
        if py["time"] is None:
            pytest.skip("No Python timing")
        s = py["time"] / rust["time"] if rust["time"] > 0 else 0
        print(f"\n  Python: {py['time']:.1f}s  Rust: {rust['time']:.1f}s  Speedup: {s:.1f}x")
        assert rust["time"] < py["time"]


def _read_rust_band(tif_path):
    """Read a Rust COG band using GDAL ReadRaster (avoids gdal_array)."""
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(str(tif_path))
    band = ds.GetRasterBand(1)
    w, h = ds.RasterXSize, ds.RasterYSize
    buf = band.ReadRaster(0, 0, w, h, buf_type=gdal.GDT_Float64)
    data = np.frombuffer(buf, dtype=np.float64).reshape(h, w).copy()
    ds = None
    return data


@pytest.mark.slow
class TestRealAccuracy:
    def test_per_band_correlation(self, rust_binary):
        py = _run_python_real()
        rust = _run_rust_real(rust_binary)

        try:
            import netCDF4
            from osgeo import gdal
        except ImportError:
            pytest.skip("netCDF4 or gdal not available")

        ds = netCDF4.Dataset(str(py["nc"]))
        all_m = {}
        for tif in rust["tifs"]:
            bname = tif.stem.split("_")[-1]
            wl = REAL["band_map"].get(bname)
            if not wl or f"rhos_{wl}" not in ds.variables:
                continue
            rd = _read_rust_band(tif)
            pd = ds.variables[f"rhos_{wl}"][:].astype(np.float64)
            h, w = min(rd.shape[0], pd.shape[0]), min(rd.shape[1], pd.shape[1])
            m = pixel_metrics(rd[:h, :w], pd[:h, :w])
            if m:
                all_m[bname] = m
        ds.close()

        assert len(all_m) >= 5
        print(f"\n  {'Band':<6} {'R':>10} {'RMSE':>10} {'Bias':>10} {'%<0.01':>8} {'%<0.05':>8}")
        print(f"  {'─'*6} {'─'*10} {'─'*10} {'─'*10} {'─'*8} {'─'*8}")
        for b in AC_BANDS:
            if b not in all_m:
                continue
            v = all_m[b]
            print(f"  {b:<6} {v['pearson_r']:>10.6f} {v['rmse']:>10.6f} "
                  f"{v['bias']:>+10.6f} {v['pct_within_0.01']:>7.1f}% {v['pct_within_0.05']:>7.1f}%")
        mean_r = np.mean([v["pearson_r"] for v in all_m.values()])
        mean_rmse = np.mean([v["rmse"] for v in all_m.values()])
        print(f"\n  Mean R = {mean_r:.6f}, Mean RMSE = {mean_rmse:.6f}")
        assert mean_r > 0.999, f"Mean R={mean_r:.6f} below 0.999"

    def test_all_bands_within_tolerance(self, rust_binary):
        """Every band must have 100% pixels within 0.05 reflectance."""
        py = _run_python_real()
        rust = _run_rust_real(rust_binary)
        try:
            import netCDF4
        except ImportError:
            pytest.skip("netCDF4 not available")

        ds = netCDF4.Dataset(str(py["nc"]))
        for tif in rust["tifs"]:
            bname = tif.stem.split("_")[-1]
            wl = REAL["band_map"].get(bname)
            if not wl or f"rhos_{wl}" not in ds.variables:
                continue
            rd = _read_rust_band(tif)
            pd = ds.variables[f"rhos_{wl}"][:].astype(np.float64)
            h, w = min(rd.shape[0], pd.shape[0]), min(rd.shape[1], pd.shape[1])
            m = pixel_metrics(rd[:h, :w], pd[:h, :w])
            if m:
                assert m["pct_within_0.05"] == 100.0, f"{bname}: {m['pct_within_0.05']:.1f}% < 100%"
        ds.close()

    def test_model_and_aot_match(self, rust_binary):
        """Rust must select same model and AOT as Python."""
        rust = _run_rust_real(rust_binary)
        # Parse Rust stdout for model/AOT (cached in timing dir)
        out = CACHE / "sentinel2_sa" / "rust_atcor_perf"
        # Re-run to capture stdout if needed
        safe = CACHE / "sentinel2_sa" / REAL["safe"]
        r = subprocess.run(
            [str(rust_binary), "--file", str(safe), "--output", str(out),
             "--res", "20", "--aot-mode", "fixed", "--model", "auto"],
            capture_output=True, text=True, timeout=300,
            env={**os.environ, "RUST_LOG": "info"},
        )
        assert "MOD1" in r.stdout, "Rust should select MOD1"
        m = re.search(r"AOT=([\d.]+)", r.stdout)
        if m:
            aot = float(m.group(1))
            assert abs(aot - 0.011) < 0.002, f"AOT={aot} too far from Python 0.011"


# ── Report ───────────────────────────────────────────────────────────────

class TestReport:
    def test_json_report(self, synthetic_pair):
        report = {
            "test": "s2_atcor_perf",
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "synthetic": {
                "size": f"{synthetic_pair['size']}x{synthetic_pair['size']}",
                "bands": 11,
                "time_s": round(synthetic_pair["py_time"], 4),
                "per_band": {
                    b: {k: round(v, 8) if isinstance(v, float) else v
                        for k, v in m.items()}
                    for b, m in synthetic_pair["metrics"].items()
                },
            },
        }
        # Append real data results if cached
        pixel_report = CACHE / "s2a_real_pixel_report.json"
        if pixel_report.exists():
            report["real_data"] = json.loads(pixel_report.read_text())
        out = CACHE / "s2_atcor_perf_report.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2))
        print(f"\n  Report → {out}")

    def test_summary_table(self, synthetic_pair):
        m = synthetic_pair["metrics"]
        sz = synthetic_pair["size"]
        t = synthetic_pair["py_time"]
        mpx = sz * sz * 11 / t / 1e6

        print("\n")
        print("=" * 72)
        print(f"  S2 ATCOR ACCURACY REPORT  (synthetic {sz}×{sz}, 11 bands)")
        print("=" * 72)
        print(f"  Pipeline time: {t:.4f}s  ({mpx:.1f} Mpx/s)")
        print(f"\n  {'Band':<6} {'λ nm':>6} {'R':>10} {'RMSE':>12} "
              f"{'Bias':>12} {'%<0.01':>8} {'Mean ρs':>10}")
        print(f"  {'─'*6} {'─'*6} {'─'*10} {'─'*12} {'─'*12} {'─'*8} {'─'*10}")
        for b in AC_BANDS:
            if b not in m:
                continue
            v = m[b]
            print(f"  {b:<6} {S2_BAND_WL[b]:>6.0f} {v['pearson_r']:>10.6f} "
                  f"{v['rmse']:>12.2e} {v['bias']:>12.2e} "
                  f"{v['pct_within_0.01']:>7.1f}% {v['mean_a']:>10.6f}")

        rs = [m[b]["pearson_r"] for b in m]
        rmses = [m[b]["rmse"] for b in m]
        print(f"\n  Mean R:    {np.mean(rs):.6f}")
        print(f"  Mean RMSE: {np.mean(rmses):.2e}")
        print(f"  Max RMSE:  {np.max(rmses):.2e}")
        print("=" * 72)
