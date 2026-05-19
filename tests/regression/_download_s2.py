#!/usr/bin/env python3
"""Download Sentinel-2 L1C SAFE from AWS S3 (requester-pays).

Reconstructs a minimal .SAFE directory with JP2 bands + XML metadata
from the sentinel-s2-l1c bucket.

Usage: python _download_s2.py S2A   # or S2B
"""

import json, os, subprocess, sys
from pathlib import Path

CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache")) / "sentinel2_sa"

SCENES = {
    "S2A": {
        "safe": "S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450.SAFE",
        "tile_path": "tiles/54/H/TF/2024/3/21/0",
        "product_path": "products/2024/3/21/S2A_MSIL1C_20240321T003701_N0510_R059_T54HTF_20240321T020450",
    },
    "S2B": {
        "safe": "S2B_MSIL1C_20240228T004659_N0510_R102_T54HTF_20240228T020941.SAFE",
        "tile_path": "tiles/54/H/TF/2024/2/28/0",
        "product_path": "products/2024/2/28/S2B_MSIL1C_20240228T004659_N0510_R102_T54HTF_20240228T020941",
    },
}

BUCKET = "sentinel-s2-l1c"
BANDS = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B09", "B10", "B11", "B12"]


AWS_PROFILE = os.environ.get("AWS_PROFILE", "adias-prod-sso-power")
AWS_REGION = "eu-central-1"


def s3_cp(s3_path, local_path):
    """Download from requester-pays S3 bucket."""
    local_path.parent.mkdir(parents=True, exist_ok=True)
    if local_path.exists() and local_path.stat().st_size > 0:
        return True
    cmd = ["aws", "s3", "cp", f"s3://{BUCKET}/{s3_path}", str(local_path),
           "--request-payer", "requester", "--region", AWS_REGION,
           "--profile", AWS_PROFILE]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        print(f"  WARN: {r.stderr.strip()}", file=sys.stderr)
        return False
    return True


def download(label):
    cfg = SCENES[label]
    safe_dir = CACHE / cfg["safe"]
    tile = cfg["tile_path"]

    # Get tileInfo.json first to discover granule structure
    info_local = safe_dir / "tileInfo.json"
    print(f"Downloading {label} tile info...")
    s3_cp(f"{tile}/tileInfo.json", info_local)

    if not info_local.exists():
        print("Failed to get tileInfo.json", file=sys.stderr)
        return False

    info = json.loads(info_local.read_text())
    product_name = info.get("productName", cfg["safe"].replace(".SAFE", ""))

    # Build SAFE directory structure
    # .SAFE/GRANULE/<granule_id>/IMG_DATA/
    # .SAFE/MTD_MSIL1C.xml
    granule_id = f"L1C_{info['path'].replace('tiles/','T').replace('/','_')}"
    # Typical: L1C_T54HTF_A045123_20240321T004321
    # Use a simpler name from tileInfo
    granule_id = f"L1C_T54HTF_{label}"

    img_dir = safe_dir / "GRANULE" / granule_id / "IMG_DATA"
    img_dir.mkdir(parents=True, exist_ok=True)
    granule_dir = safe_dir / "GRANULE" / granule_id

    # Download bands
    print(f"Downloading {len(BANDS)} bands...")
    for band in BANDS:
        local = img_dir / f"T54HTF_20240321T003701_{band}.jp2"
        if not s3_cp(f"{tile}/{band}.jp2", local):
            print(f"  Failed: {band}")

    # Download granule metadata (MTD_TL.xml)
    print("Downloading metadata...")
    s3_cp(f"{tile}/metadata.xml", granule_dir / "MTD_TL.xml")

    # Download product-level metadata
    # Try to get MTD_MSIL1C.xml from the product path
    prod = cfg["product_path"]
    s3_cp(f"{prod}/metadata.xml", safe_dir / "MTD_MSIL1C.xml")

    # If product metadata failed, create a minimal one
    mtd = safe_dir / "MTD_MSIL1C.xml"
    if not mtd.exists() or mtd.stat().st_size == 0:
        print("Creating minimal MTD_MSIL1C.xml...")
        mtd.write_text(f"""<?xml version="1.0" encoding="UTF-8"?>
<n1:Level-1C_User_Product>
  <n1:General_Info>
    <Product_Info>
      <PRODUCT_START_TIME>{info.get('timestamp','2024-03-21T00:47:04Z')}</PRODUCT_START_TIME>
      <SPACECRAFT_NAME>Sentinel-2A</SPACECRAFT_NAME>
      <PROCESSING_BASELINE>05.10</PROCESSING_BASELINE>
    </Product_Info>
    <Product_Image_Characteristics>
      <QUANTIFICATION_VALUE>10000</QUANTIFICATION_VALUE>
      <Radiometric_Offset_List>
        <RADIO_ADD_OFFSET band_id="0">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="1">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="2">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="3">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="4">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="5">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="6">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="7">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="8">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="9">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="10">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="11">-1000</RADIO_ADD_OFFSET>
        <RADIO_ADD_OFFSET band_id="12">-1000</RADIO_ADD_OFFSET>
      </Radiometric_Offset_List>
    </Product_Image_Characteristics>
  </n1:General_Info>
</n1:Level-1C_User_Product>""")

    # Verify
    jp2s = list(img_dir.glob("*.jp2"))
    print(f"\nDone: {len(jp2s)} JP2 files in {safe_dir}")
    return len(jp2s) >= 13


if __name__ == "__main__":
    label = sys.argv[1] if len(sys.argv) > 1 else "S2A"
    if label not in SCENES:
        print(f"Unknown label: {label}. Use S2A or S2B.")
        sys.exit(1)
    ok = download(label)
    sys.exit(0 if ok else 1)
