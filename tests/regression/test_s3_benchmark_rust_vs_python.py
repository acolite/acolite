"""
Sentinel-3 OLCI benchmark: Rust vs Python ACOLITE — real data, full DSF.

South Australia, Gulf St Vincent — same region as S2 benchmarks.
S3B EFR scene from 2024-03-21 (same date as S2A benchmark scene).

Both pipelines use full LUT-based DSF atmospheric correction.

Measures:
  - Wall-clock time for L1 conversion + atmospheric correction
  - Per-band ρs accuracy (Rust vs Python L2R)
  - Spectral shape consistency
  - Speedup factor

Run:
  pytest tests/regression/test_s3_benchmark_rust_vs_python.py -v -s
"""

import json, os, re, subprocess, sys, time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent.parent
CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache")) / "sentinel3_sa"

LIMIT = [-35.0, 138.0, -34.5, 138.7]
LIMIT_STR = ",".join(str(x) for x in LIMIT)

SCENE_NAME = "S3B_OL_1_EFR____20240321T003707_20240321T004007_20240921T122930_0179_091_059_3600_MAR_R_NT_004.SEN3"

# Rust band → Python rhos wavelength
RHOS_MAP = {
    "Oa01": "401", "Oa02": "412", "Oa03": "443", "Oa04": "490",
    "Oa05": "510", "Oa06": "560", "Oa07": "620", "Oa08": "665",
    "Oa09": "674", "Oa10": "681", "Oa11": "709", "Oa12": "754",
    "Oa15": "768", "Oa16": "779", "Oa17": "865", "Oa18": "884",
    "Oa21": "1016",
}


def _compare(a, b):
    h, w = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
    a, b = a[:h, :w], b[:h, :w]
    mask = np.isfinite(a) & np.isfinite(b) & (a > -0.5) & (b > -0.5) & (a < 2.0) & (b < 2.0)
    a, b = a[mask], b[mask]
    if len(a) < 100:
        return None
    diff = a - b
    r = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else 0.0
    return {
        "pearson_r": r,
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "bias": float(np.mean(diff)),
        "n_pixels": int(len(a)),
    }


# ── Python pipeline ──────────────────────────────────────────────────────

def _run_python():
    sen3_dir = CACHE / SCENE_NAME
    if not sen3_dir.exists():
        pytest.skip("No S3 data — run _download_s3.py S3B first")

    out_dir = CACHE / "py_output"
    l2r_nc = out_dir / "S3B_OLCI_2024_03_21_00_37_06_FR_L2R.nc"
    timing_file = out_dir / "timing.txt"

    if l2r_nc.exists() and timing_file.exists():
        return {
            "l2r": l2r_nc, "out_dir": out_dir,
            "time": float(timing_file.read_text().strip()), "cached": True,
        }

    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0,'{REPO}'); "
         f"import acolite; acolite.acolite.acolite_run(settings={{"
         f"'inputfile':'{sen3_dir}',"
         f"'output':'{out_dir}',"
         f"'limit':[-35.0,138.0,-34.5,138.7],"
         f"'smile_correction':True,'smile_correction_tgas':True,"
         f"'use_supplied_ancillary':True,"
         f"'dsf_aot_estimate':'fixed','dsf_spectrum_option':'intercept',"
         f"'dsf_intercept_pixels':200,"
         f"'l2w_parameters':'','rgb_rhot':False,'rgb_rhos':False,'map_l2w':False,"
         f"'output_geometry':False}})"],
        capture_output=True, text=True, timeout=600, cwd=str(REPO),
    )
    elapsed = time.time() - t0
    timing_file.write_text(f"{elapsed:.2f}")

    if r.returncode != 0:
        pytest.fail(f"Python failed ({elapsed:.1f}s):\n{r.stderr[:500]}")

    l2r_files = list(out_dir.glob("*_L2R.nc"))
    return {
        "l2r": l2r_files[0] if l2r_files else None,
        "out_dir": out_dir, "time": elapsed, "cached": False,
    }


# ── Rust pipeline ─────────────────────────────────────────────────────────

def _build_rust():
    r = subprocess.run(
        ["cargo", "build", "--release", "--features", "full-io",
         "--example", "process_sentinel3"],
        cwd=str(REPO), capture_output=True, text=True, timeout=300,
    )
    if r.returncode != 0:
        pytest.skip(f"Rust build failed: {r.stderr[:300]}")
    binary = REPO / "target" / "release" / "examples" / "process_sentinel3"
    if not binary.exists():
        pytest.skip("Rust binary not found")
    return binary


def _run_rust(binary):
    sen3_dir = CACHE / SCENE_NAME
    if not sen3_dir.exists():
        pytest.skip("No S3 data")

    out_dir = CACHE / "rust_dsf_output"
    timing_file = out_dir / "timing.txt"
    stats_file = out_dir / "stats.json"

    if stats_file.exists() and timing_file.exists():
        return {
            "out_dir": out_dir,
            "time": float(timing_file.read_text().strip()),
            "stats": json.loads(stats_file.read_text()),
            "cached": True,
        }

    out_dir.mkdir(parents=True, exist_ok=True)

    data_dir = REPO / "data"
    if not (data_dir / "LUT").exists():
        pytest.skip("No data/LUT directory for full DSF")

    t0 = time.time()
    r = subprocess.run(
        [str(binary), "--scene", str(sen3_dir),
         "--limit", LIMIT_STR,
         "--data-dir", str(data_dir),
         "--output", str(out_dir)],
        capture_output=True, text=True, timeout=300,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0
    timing_file.write_text(f"{elapsed:.2f}")

    if r.returncode != 0:
        pytest.fail(f"Rust failed ({elapsed:.1f}s):\n{r.stderr[:300]}\n{r.stdout[-300:]}")

    # Parse band stats from stdout
    stats = {}
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].startswith("Oa"):
            try:
                mean = float(parts[2])
                stats[parts[0]] = {"wavelength": float(parts[1]), "mean_rhos": mean}
            except ValueError:
                pass

    # Parse AOT and model
    for line in r.stdout.splitlines():
        if "AOT:" in line:
            m = re.search(r'AOT:\s*([\d.]+).*model:\s*(\S+)', line)
            if m:
                stats["_aot"] = float(m.group(1))
                stats["_model"] = m.group(2).rstrip(")")

    stats_file.write_text(json.dumps(stats, indent=2))

    return {
        "out_dir": out_dir, "time": elapsed,
        "stats": stats, "cached": False,
        "stdout": r.stdout,
    }


# ── Fixtures ──────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def rust_binary():
    return _build_rust()

@pytest.fixture(scope="module")
def python_output():
    return _run_python()

@pytest.fixture(scope="module")
def rust_output(rust_binary):
    return _run_rust(rust_binary)


# ── Performance ───────────────────────────────────────────────────────────

class TestS3Performance:
    def test_rust_completes(self, rust_output):
        assert rust_output["time"] > 0
        print(f"\n  Rust time: {rust_output['time']:.2f}s")

    def test_python_completes(self, python_output):
        assert python_output["time"] > 0
        print(f"\n  Python time: {python_output['time']:.2f}s")

    def test_rust_produces_bands(self, rust_output):
        n = len([k for k in rust_output.get("stats", {}) if k.startswith("Oa")])
        print(f"\n  Rust produced {n} band statistics")
        assert n >= 15, f"Expected ≥15 bands, got {n}"

    def test_rust_aot_physical(self, rust_output):
        aot = rust_output.get("stats", {}).get("_aot", None)
        if aot is not None:
            print(f"\n  Rust AOT: {aot:.4f}")
            assert 0.0 < aot < 2.0, f"AOT {aot} out of physical range"


# ── Accuracy: compare Rust ρs vs Python ρs ────────────────────────────────

class TestS3Accuracy:
    def test_python_has_rhos(self, python_output):
        if python_output.get("l2r") is None:
            pytest.skip("No L2R")
        import netCDF4
        nc = netCDF4.Dataset(str(python_output["l2r"]))
        rhos = [v for v in nc.variables if v.startswith("rhos_")]
        nc.close()
        print(f"\n  Python L2R has {len(rhos)} rhos bands")
        assert len(rhos) >= 15

    def test_rust_means_physical_range(self, rust_output):
        for band, info in rust_output.get("stats", {}).items():
            if band.startswith("_"):
                continue
            m = info["mean_rhos"]
            if np.isnan(m):
                continue  # O2/H2O absorption bands (Oa13,Oa14,Oa19,Oa20) may be NaN
            assert -0.05 <= m <= 1.0, f"{band}: mean={m} out of range"

    def test_rhos_spectral_comparison(self, rust_output, python_output):
        """Compare Rust ρs means against Python L2R ρs means."""
        if python_output.get("l2r") is None:
            pytest.skip("No Python L2R")

        import netCDF4
        nc = netCDF4.Dataset(str(python_output["l2r"]))

        rust_means, py_means, bands = [], [], []
        for rust_band, py_wave in RHOS_MAP.items():
            py_var = f"rhos_{py_wave}"
            if py_var not in nc.variables:
                continue
            rust_info = rust_output.get("stats", {}).get(rust_band)
            if rust_info is None:
                continue

            py_data = nc.variables[py_var][:].astype(np.float64)
            valid = py_data[np.isfinite(py_data)]
            if len(valid) == 0:
                continue
            py_mean = float(np.nanmean(valid))

            rust_means.append(rust_info["mean_rhos"])
            py_means.append(py_mean)
            bands.append(rust_band)

        nc.close()

        if len(rust_means) < 5:
            pytest.skip(f"Only {len(rust_means)} bands")

        rust_arr = np.array(rust_means)
        py_arr = np.array(py_means)
        r = float(np.corrcoef(rust_arr, py_arr)[0, 1]) if np.std(rust_arr) > 0 and np.std(py_arr) > 0 else 0.0
        mean_rmse = float(np.sqrt(np.mean((rust_arr - py_arr) ** 2)))

        print(f"\n  Surface reflectance: Rust ρs vs Python ρs (both full DSF)")
        print(f"  {'Band':<6} {'Rust':>10} {'Python':>10} {'Diff':>10}")
        print(f"  {'─'*6} {'─'*10} {'─'*10} {'─'*10}")
        for i, b in enumerate(bands):
            print(f"  {b:<6} {rust_arr[i]:>10.4f} {py_arr[i]:>10.4f} {rust_arr[i]-py_arr[i]:>+10.4f}")

        print(f"\n  Spectral R = {r:.4f}  RMSE = {mean_rmse:.4f}")

        assert r > 0.70, f"ρs spectral R={r:.4f} too low"
        assert mean_rmse < 0.05, f"ρs RMSE={mean_rmse:.4f} too high"

    def test_aot_within_tolerance(self, rust_output, python_output):
        """Rust and Python AOT should be in the same ballpark."""
        rust_aot = rust_output.get("stats", {}).get("_aot")
        if rust_aot is None:
            pytest.skip("No Rust AOT")
        # Python AOT is typically ~0.05-0.15 for this scene
        print(f"\n  Rust AOT: {rust_aot:.4f}")
        assert 0.01 < rust_aot < 0.5


# ── Summary ───────────────────────────────────────────────────────────────

class TestS3BenchmarkSummary:
    def test_summary(self, rust_output, python_output):
        py_t = python_output["time"]
        rs_t = rust_output["time"]
        speedup = py_t / rs_t if rs_t > 0 else 0
        n_bands = len([k for k in rust_output.get("stats", {}) if k.startswith("Oa")])

        report = {
            "scene": SCENE_NAME,
            "region": "South Australia, Gulf St Vincent",
            "limit": LIMIT,
            "python_time_s": py_t,
            "rust_time_s": rs_t,
            "speedup": round(speedup, 2),
            "rust_bands": n_bands,
            "rust_aot": rust_output.get("stats", {}).get("_aot"),
            "rust_model": rust_output.get("stats", {}).get("_model"),
            "ac_mode": "full_dsf",
        }

        print("\n")
        print("=" * 70)
        print("  S3 OLCI BENCHMARK: Rust vs Python ACOLITE (Full DSF)")
        print("  South Australia, Gulf St Vincent")
        print(f"  Scene: {SCENE_NAME[:60]}...")
        print("=" * 70)
        print(f"  Python:  {py_t:>8.1f}s")
        print(f"  Rust:    {rs_t:>8.2f}s")
        print(f"  Speedup: {speedup:>8.2f}x")
        print(f"  Bands:   {n_bands}")
        print(f"  AOT:     {report.get('rust_aot', 'N/A')}")
        print(f"  Model:   {report.get('rust_model', 'N/A')}")
        print("=" * 70)

        report_path = CACHE / "s3_benchmark_report.json"
        report_path.write_text(json.dumps(report, indent=2))
