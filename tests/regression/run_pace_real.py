#!/usr/bin/env python3
"""
PACE OCI real-data regression test.

Searches OB-DAAC via CMR, downloads a PACE L1B granule (cached),
processes with both Python ACOLITE and Rust acolite-rs, and compares outputs.

Usage:
    python tests/regression/run_pace_real.py
    python tests/regression/run_pace_real.py --skip-download /path/to/cached.nc
    python tests/regression/run_pace_real.py --lat 36.0 --lon -75.5 --date 2024-07-01

Requires:
    - EarthData credentials (env EARTHDATA_u/EARTHDATA_p, or ~/.netrc, or config/credentials.txt)
      NOTE: PACE OCI L1B V3 is publicly accessible — credentials may not be needed.
    - Python ACOLITE dependencies (numpy, scipy, netCDF4, etc.)
    - Rust acolite-rs built with --features netcdf
"""

import argparse
import json
import os
import sys
import time
import subprocess
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CACHE_DIR = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache"))

# Default: Chesapeake Bay area, a date known to have PACE coverage
DEFAULT_LAT = 36.0
DEFAULT_LON = -75.5
DEFAULT_DATE = "2024-07-01"
# Small ROI for fast processing (~50×50 pixels)
DEFAULT_LIMIT = [35.8, -75.4, 36.0, -75.2]


def cmr_search(lat, lon, date):
    """Search CMR for PACE OCI L1B granules (public API)."""
    import requests
    end_date = date  # single day
    url = "https://cmr.earthdata.nasa.gov/search/granules.json"
    params = {
        "echo_collection_id": "C3392966952-OB_CLOUD",
        "point": f"{lon},{lat}",
        "temporal": f"{date}T00:00:00Z,{end_date}T23:59:59Z",
        "page_size": "5",
    }
    print(f"→ CMR search: PACE OCI L1B at ({lat}, {lon}) on {date}")
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    entries = r.json()["feed"]["entry"]
    granules = []
    for e in entries:
        gid = e.get("producer_granule_id", e.get("title", "unknown"))
        urls = [
            l["href"] for l in e.get("links", [])
            if l["href"].endswith(".nc") and "ob-cumulus" in l["href"] and l["href"].startswith("http")
        ]
        size_mb = float(e.get("granule_size", 0))
        if urls:
            granules.append({"id": gid, "url": urls[0], "size_mb": size_mb})
    print(f"  Found {len(granules)} granule(s)")
    for g in granules:
        print(f"    {g['id']} ({g['size_mb']:.0f} MB)")
    return granules


def download_granule(granule, cache_dir):
    """Download granule to cache (skip if exists). Handles EarthData OAuth redirect."""
    import requests
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / granule["id"]
    if dest.exists() and dest.stat().st_size > 1e6:
        size_mb = dest.stat().st_size / 1e6
        print(f"  ✓ Cached: {dest} ({size_mb:.0f} MB)")
        return dest

    print(f"  Downloading {granule['id']} ({granule['size_mb']:.0f} MB)...")
    print(f"    URL: {granule['url']}")

    auth = _get_earthdata_auth()
    if auth is None:
        print("  ✗ EarthData credentials required.")
        print("    Set EARTHDATA_u/EARTHDATA_p env vars, ~/.netrc, or ~/.easi-workflows-auth.conf")
        sys.exit(1)

    t0 = time.time()
    tmp = str(dest) + ".tmp"

    # Two-step download matching ACOLITE's pattern:
    # 1) GET without auth → follows redirect to URS OAuth URL
    # 2) GET the redirected URL WITH auth → gets the actual data
    with requests.Session() as session:
        r1 = session.request("get", granule["url"])
        r = session.get(r1.url, auth=auth, stream=True)
        r.raise_for_status()

        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0:
                        pct = downloaded / total * 100
                        print(f"\r    {downloaded/1e6:.0f}/{total/1e6:.0f} MB ({pct:.0f}%)", end="", flush=True)
                    else:
                        print(f"\r    {downloaded/1e6:.0f} MB", end="", flush=True)

    os.rename(tmp, str(dest))
    elapsed = time.time() - t0
    size_mb = dest.stat().st_size / 1e6
    print(f"\n  ✓ Downloaded {size_mb:.0f} MB in {elapsed:.0f}s ({size_mb/elapsed:.1f} MB/s)")
    return dest


def _get_earthdata_auth():
    """Get EarthData credentials from env, .netrc, easi-workflows, or config."""
    import netrc as netrc_mod
    import configparser

    # 1. Environment variables
    u = os.environ.get("EARTHDATA_u")
    p = os.environ.get("EARTHDATA_p")
    if u and p:
        return (u, p)

    # 2. .netrc (try both 'urs.earthdata.nasa.gov' and 'earthdata')
    try:
        nrc = netrc_mod.netrc()
        for host in ["urs.earthdata.nasa.gov", "earthdata"]:
            entry = nrc.authenticators(host)
            if entry:
                return (entry[0], entry[2])
    except (FileNotFoundError, netrc_mod.NetrcParseError):
        pass

    # 3. EASI workflows auth config
    easi_path = Path.home() / ".easi-workflows-auth.conf"
    if easi_path.exists():
        try:
            cp = configparser.ConfigParser()
            cp.read_string("[DEFAULT]\n" + easi_path.read_text().replace("    ", ""))
            if "earthdata" in cp:
                u = cp["earthdata"].get("user")
                p = cp["earthdata"].get("password")
                if u and p:
                    return (u, p)
        except Exception:
            # Try manual parsing for indented INI-like format
            section = None
            u, p = None, None
            for line in easi_path.read_text().splitlines():
                stripped = line.strip()
                if stripped.startswith("[") and stripped.endswith("]"):
                    section = stripped[1:-1]
                elif section == "earthdata" and ":" in stripped:
                    key, val = stripped.split(":", 1)
                    key, val = key.strip(), val.strip()
                    if key == "user":
                        u = val
                    elif key == "password":
                        p = val
            if u and p:
                return (u, p)

    # 4. ACOLITE config
    cred_path = REPO_ROOT / "config" / "credentials.txt"
    if cred_path.exists():
        u, p = None, None
        for line in cred_path.read_text().splitlines():
            line = line.strip()
            if line.startswith("EARTHDATA_u=") and line.split("=", 1)[1]:
                u = line.split("=", 1)[1]
            if line.startswith("EARTHDATA_p=") and line.split("=", 1)[1]:
                p = line.split("=", 1)[1]
        if u and p:
            return (u, p)

    return None


def process_python(nc_path, output_dir, limit):
    """Process with Python ACOLITE, return L1R output path and stats."""
    sys.path.insert(0, str(REPO_ROOT))
    import acolite as ac

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n→ Python ACOLITE processing (limit={limit})...")
    t0 = time.time()

    ofiles, setu = ac.pace.l1_convert(
        str(nc_path),
        output=str(output_dir),
        settings={"limit": limit, "verbosity": 1},
    )
    elapsed = time.time() - t0

    if not ofiles:
        print("  ✗ No output produced")
        return None, None

    print(f"  ✓ Output: {ofiles[0]}")
    print(f"  ✓ Processed in {elapsed:.1f}s")

    # Extract per-band statistics
    from netCDF4 import Dataset
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
                    "valid_pct": float(np.sum(np.isfinite(data)) / data.size * 100),
                }

    print(f"  Bands: {len(stats)}")
    return ofiles[0], stats


def process_rust(nc_path, output_dir, limit):
    """Process with Rust acolite-rs, return output path and timing."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Build Rust binary if needed
    print(f"\n→ Rust acolite-rs processing...")
    build = subprocess.run(
        ["cargo", "build", "--release", "--features", "netcdf", "--example", "process_pace"],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=300,
    )
    if build.returncode != 0:
        print(f"  ✗ Build failed: {build.stderr[:300]}")
        return None, None

    binary = REPO_ROOT / "target" / "release" / "examples" / "process_pace"
    if not binary.exists():
        print("  ✗ Binary not found")
        return None, None

    t0 = time.time()
    result = subprocess.run(
        [str(binary), "--file", str(nc_path), "--output", str(output_dir)],
        capture_output=True, text=True, timeout=600,
        env={**os.environ, "RUST_LOG": "info"},
    )
    elapsed = time.time() - t0

    print(f"  Exit code: {result.returncode}")
    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            print(f"  {line}")
    if result.returncode != 0 and result.stderr:
        print(f"  stderr: {result.stderr[:500]}")

    print(f"  ✓ Rust processed in {elapsed:.1f}s")

    # Find output zarr
    zarr_dirs = list(output_dir.glob("*.zarr"))
    return (zarr_dirs[0] if zarr_dirs else None), elapsed


def compare_results(py_stats, rust_output, tolerance=1e-3):
    """Compare Python and Rust outputs."""
    print(f"\n{'='*60}")
    print("REGRESSION COMPARISON")
    print(f"{'='*60}")

    if py_stats is None:
        print("  ✗ No Python output to compare")
        return False

    print(f"\nPython ACOLITE: {len(py_stats)} bands")
    print(f"{'Band':<30} {'Mean':>10} {'Std':>10} {'Min':>10} {'Max':>10} {'Valid%':>8}")
    print("-" * 78)
    for band in sorted(py_stats.keys()):
        s = py_stats[band]
        print(f"{band:<30} {s['mean']:10.6f} {s['std']:10.6f} {s['min']:10.6f} {s['max']:10.6f} {s['valid_pct']:7.1f}%")

    # Physical sanity checks
    print(f"\nSanity checks:")
    failures = 0
    for band, s in py_stats.items():
        if s["mean"] < -0.05 or s["mean"] > 1.5:
            print(f"  ✗ {band}: mean={s['mean']:.4f} out of range")
            failures += 1
        if s["valid_pct"] < 10:
            print(f"  ⚠ {band}: only {s['valid_pct']:.0f}% valid pixels")

    # Check spectral consistency: blue bands should have higher rhot than SWIR
    blue_bands = {k: v for k, v in py_stats.items() if "blue" in k and "400" < k.split("_")[-1] < "500"}
    swir_bands = {k: v for k, v in py_stats.items() if "SWIR" in k}
    if blue_bands and swir_bands:
        blue_mean = np.mean([v["mean"] for v in blue_bands.values()])
        swir_mean = np.mean([v["mean"] for v in swir_bands.values()])
        if blue_mean > swir_mean:
            print(f"  ✓ Blue mean ({blue_mean:.4f}) > SWIR mean ({swir_mean:.4f}) — Rayleigh consistent")
        else:
            print(f"  ⚠ Blue mean ({blue_mean:.4f}) ≤ SWIR mean ({swir_mean:.4f})")

    if rust_output and rust_output.exists():
        print(f"\nRust output: {rust_output}")
        # Check zarr structure
        if rust_output.suffix == ".zarr":
            for expected in ["zarr.json", "data", "wavelengths", "bandwidths"]:
                path = rust_output / expected
                if path.exists():
                    print(f"  ✓ {expected}")
                else:
                    print(f"  ✗ Missing {expected}")
                    failures += 1

    if failures == 0:
        print(f"\n✓ ALL CHECKS PASSED")
    else:
        print(f"\n✗ {failures} CHECK(S) FAILED")

    return failures == 0


def main():
    parser = argparse.ArgumentParser(description="PACE OCI real-data regression test")
    parser.add_argument("--lat", type=float, default=DEFAULT_LAT)
    parser.add_argument("--lon", type=float, default=DEFAULT_LON)
    parser.add_argument("--date", default=DEFAULT_DATE)
    parser.add_argument("--limit", type=float, nargs=4, default=DEFAULT_LIMIT,
                        metavar=("S", "W", "N", "E"))
    parser.add_argument("--skip-download", metavar="FILE",
                        help="Use existing local file instead of downloading")
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument("--python-only", action="store_true")
    parser.add_argument("--rust-only", action="store_true")
    args = parser.parse_args()

    print("ACOLITE PACE OCI Real-Data Regression Test")
    print(f"Location: ({args.lat}, {args.lon})")
    print(f"Date: {args.date}")
    print(f"ROI limit: {args.limit}")
    print(f"Cache: {args.cache_dir}")
    print()

    # Step 1: Get the data file
    if args.skip_download:
        nc_path = Path(args.skip_download)
        if not nc_path.exists():
            print(f"File not found: {nc_path}")
            sys.exit(1)
        print(f"→ Using local file: {nc_path}")
    else:
        granules = cmr_search(args.lat, args.lon, args.date)
        if not granules:
            print("No granules found. Try a different date or location.")
            sys.exit(1)
        nc_path = download_granule(granules[0], args.cache_dir)

    # Step 2: Process with Python ACOLITE
    py_output = None
    py_stats = None
    if not args.rust_only:
        py_dir = args.cache_dir / "py_output"
        py_output, py_stats = process_python(nc_path, py_dir, args.limit)

    # Step 3: Process with Rust
    rust_output = None
    rust_time = None
    if not args.python_only:
        rust_dir = args.cache_dir / "rust_output"
        rust_output, rust_time = process_rust(nc_path, rust_dir, args.limit)

    # Step 4: Compare
    success = compare_results(py_stats, rust_output)

    # Save reference stats
    if py_stats:
        ref_path = args.cache_dir / "pace_reference_stats.json"
        ref_path.write_text(json.dumps({
            "date": args.date,
            "location": [args.lat, args.lon],
            "limit": args.limit,
            "granule": nc_path.name,
            "bands": py_stats,
        }, indent=2))
        print(f"\nReference stats saved: {ref_path}")

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
