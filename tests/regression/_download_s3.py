#!/usr/bin/env python3
"""Download Sentinel-3 OLCI EFR from Copernicus Data Space (CDSE).

Requires CDSE credentials in ~/.netrc (machine cdse) or environment.
See: https://dataspace.copernicus.eu/

Usage: python _download_s3.py S3B
"""

import os, subprocess, sys, zipfile
from pathlib import Path

CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache")) / "sentinel3_sa"
REPO = Path(__file__).resolve().parent.parent.parent

SCENES = {
    "S3B": {
        "scene": "S3B_OL_1_EFR____20240321T003707_20240321T004007_20240921T122930_0179_091_059_3600_MAR_R_NT_004.SEN3",
    },
    "S3A": {
        "scene": "S3A_OL_1_EFR____20240320T000118_20240320T000418_20240920T113545_0179_110_187_3600_MAR_R_NT_004.SEN3",
    },
}


def download(label):
    cfg = SCENES[label]
    scene_name = cfg["scene"]
    sen3_dir = CACHE / scene_name

    # Already downloaded?
    if sen3_dir.exists():
        ncs = list(sen3_dir.glob("*.nc"))
        if len(ncs) >= 22:  # 21 radiance + geo + instrument + ...
            print(f"Already cached: {sen3_dir} ({len(ncs)} NC files)")
            return True

    CACHE.mkdir(parents=True, exist_ok=True)

    # Use ACOLITE's CDSE API
    sys.path.insert(0, str(REPO))
    import acolite as ac

    print(f"Querying CDSE for {scene_name}...")
    urls, scenes = ac.api.cdse.query(scene=scene_name, verbosity=0)
    if not urls:
        print("Scene not found on CDSE", file=sys.stderr)
        return False

    print(f"Downloading {scenes[0]}...")
    try:
        ac.api.cdse.download(urls[:1], scenes=scenes[:1], output=str(CACHE), verbosity=1)
    except Exception as e:
        print(f"Download failed: {e}", file=sys.stderr)
        return False

    # Extract if zip
    zips = list(CACHE.glob(f"*{label}*.zip"))
    for z in zips:
        print(f"Extracting {z.name}...")
        with zipfile.ZipFile(z) as zf:
            zf.extractall(CACHE)
        z.unlink()

    # Verify
    if sen3_dir.exists():
        ncs = list(sen3_dir.glob("*.nc"))
        print(f"Done: {len(ncs)} NC files in {sen3_dir}")
        return len(ncs) >= 22
    else:
        # Check if extracted with different name
        sen3_dirs = list(CACHE.glob("*.SEN3"))
        for d in sen3_dirs:
            ncs = list(d.glob("*.nc"))
            if len(ncs) >= 22:
                print(f"Found: {d} ({len(ncs)} NC files)")
                return True
        print("No valid SEN3 directory found after download", file=sys.stderr)
        return False


if __name__ == "__main__":
    label = sys.argv[1] if len(sys.argv) > 1 else "S3B"
    if label not in SCENES:
        print(f"Unknown label: {label}. Use S3A or S3B.")
        sys.exit(1)
    ok = download(label)
    sys.exit(0 if ok else 1)
