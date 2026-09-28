#!/usr/bin/env python3
"""
End-to-end integration test: Sentinel-3 OLCI → ACOLITE DSF → COG

Searches CMR for a recent S3A OLCI EFR scene over SE Australia,
downloads, processes with ACOLITE DSF, reprojects to UTM, and
exports as per-band GeoTIFFs.

Designed to run in CI (GitHub Actions) with:
  - ~15 min timeout (download + AC + export)
  - ~3 GB RAM peak
  - ~2 GB disk for scene + outputs
  - EarthData credentials via secrets

Usage:
  pytest tests/regression/test_s3_e2e_pipeline.py -v --timeout=900

Environment variables:
  EARTHDATA_u / EARTHDATA_p: EarthData credentials (or .netrc)
  ACOLITE_TEST_CACHE: Cache directory for downloaded scenes (default: /tmp/acolite_s3_cache)
"""

import os
import sys
import json
import time
import resource
import zipfile
import tempfile
from pathlib import Path

import numpy as np
import pytest
import requests

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

CACHE = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_s3_cache"))

# ─── Configuration ────────────────────────────────────────────────────────────
# SE Australia / Western Victoria coast — known to have frequent S3A coverage
# This bbox reliably intersects OLCI EFR descending passes
SEARCH_BBOX = [141.0, -36.0, 145.0, -34.0]  # [west, south, east, north] for CMR
ACOLITE_LIMIT = [-36.0, 141.0, -34.5, 143.0]  # [south, west, north, east] for ACOLITE

CMR_SEARCH_URL = "https://cmr.earthdata.nasa.gov/search/granules.json"
COLLECTION_SHORT_NAME = "S3A_OL_1_EFR"
PROVIDER = "LAADS"

# Expected outputs
MIN_COG_FILES = 20  # At minimum we expect 15 rhos + some rhot bands
MIN_COG_SIZE_KB = 100  # Each COG should be > 100 KB

# Performance thresholds
MAX_PROCESSING_TIME_S = 600  # 10 min max for AC + reprojection
MAX_PEAK_RSS_MB = 4096  # 4 GB max


# ─── Fixtures ─────────────────────────────────────────────────────────────────

def get_earthdata_credentials():
    """Get EarthData credentials from env or .netrc."""
    u = os.environ.get("EARTHDATA_u")
    p = os.environ.get("EARTHDATA_p")
    if u and p:
        return u, p

    try:
        import netrc
        nrc = netrc.netrc()
        auth = nrc.authenticators("urs.earthdata.nasa.gov") or nrc.authenticators("earthdata")
        if auth:
            return auth[0], auth[2]
    except Exception:
        pass

    pytest.skip("No EarthData credentials (set EARTHDATA_u/EARTHDATA_p or .netrc)")


def search_scene():
    """Find a recent S3A OLCI scene over SE Australia."""
    from datetime import datetime, timedelta

    end = datetime.utcnow()
    start = end - timedelta(days=14)

    params = {
        "short_name": COLLECTION_SHORT_NAME,
        "provider": PROVIDER,
        "bounding_box": f"{SEARCH_BBOX[0]},{SEARCH_BBOX[1]},{SEARCH_BBOX[2]},{SEARCH_BBOX[3]}",
        "temporal": f"{start.strftime('%Y-%m-%dT00:00:00Z')},{end.strftime('%Y-%m-%dT23:59:59Z')}",
        "page_size": "10",
        "sort_key[]": "-start_date",
    }

    resp = requests.get(CMR_SEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    entries = resp.json()["feed"]["entry"]

    for entry in entries:
        links = [
            l["href"] for l in entry.get("links", [])
            if l["href"].startswith("https://") and l["href"].endswith(".zip")
        ]
        if links:
            return {
                "id": entry.get("producer_granule_id", entry.get("title", "")),
                "url": links[0],
            }

    pytest.skip("No S3A OLCI scenes found in the last 14 days over SE Australia")


def download_scene(granule, credentials):
    """Download and extract S3 OLCI scene, using cache.
    
    LAADS DAAC uses OAuth2 redirect flow. We handle this by:
    1. Setting up a session that follows redirects
    2. Using .netrc for auth (works with redirect-based auth)
    3. Falling back to token-based auth if needed
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    sen3_name = granule["id"].replace(".zip", ".SEN3")
    sen3_dir = CACHE / sen3_name

    if sen3_dir.exists() and any(sen3_dir.glob("*.nc")):
        return sen3_dir

    zip_path = CACHE / granule["id"]

    # Setup session with proper auth for LAADS DAAC OAuth redirects
    session = requests.Session()
    username, password = credentials

    # Method 1: Try with .netrc (if configured, handles redirects automatically)
    try:
        import netrc as netrc_mod
        nrc = netrc_mod.netrc()
        # .netrc handles auth transparently with redirects
    except Exception:
        pass

    # Use auth tuple — requests handles redirect auth properly
    session.auth = (username, password)

    resp = session.get(granule["url"], allow_redirects=True, stream=True, timeout=120)

    # If 401, try getting a bearer token
    if resp.status_code == 401:
        token_resp = session.post(
            "https://urs.earthdata.nasa.gov/api/users/token",
            auth=(username, password),
        )
        if token_resp.ok:
            token = token_resp.json().get("access_token", "")
            session.headers["Authorization"] = f"Bearer {token}"
            session.auth = None
            resp = session.get(granule["url"], allow_redirects=True, stream=True, timeout=120)

    resp.raise_for_status()

    with open(zip_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=128 * 1024):
            f.write(chunk)

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(CACHE)
    zip_path.unlink()

    sen3_dirs = list(CACHE.glob("*.SEN3"))
    assert sen3_dirs, "No .SEN3 directory after extraction"
    return sen3_dirs[0]


# ─── Test ─────────────────────────────────────────────────────────────────────

class TestS3OlciEndToEnd:
    """End-to-end Sentinel-3 OLCI processing test."""

    @pytest.fixture(scope="class")
    def scene_dir(self):
        """Download scene (cached)."""
        creds = get_earthdata_credentials()
        granule = search_scene()
        print(f"\n  Scene: {granule['id']}")
        return download_scene(granule, creds)

    @pytest.fixture(scope="class")
    def processing_result(self, scene_dir):
        """Run ACOLITE DSF + reprojection."""
        import acolite as ac

        out_dir = CACHE / "output"
        out_dir.mkdir(parents=True, exist_ok=True)

        settings = {
            "inputfile": str(scene_dir),
            "output": str(out_dir),
            "limit": ACOLITE_LIMIT,
            "verbosity": 2,
            "output_projection": True,
            "output_projection_epsg": 32754,
            "output_projection_resolution": 300,
            "dsf_aot_estimate": "fixed",
            "l2w_parameters": None,
            "rgb_rhot": False,
            "rgb_rhos": False,
            "map_l2w": False,
        }

        t0 = time.time()
        result = ac.acolite.acolite_run(settings=settings)
        ac_time = time.time() - t0

        # Find L2R output
        l2r_files = []
        for key in result.values():
            if "l2r" in key:
                l2r_files.extend(key["l2r"])

        assert l2r_files, "No L2R output produced"
        l2r_nc = l2r_files[0]

        # Reproject
        t1 = time.time()
        projected_nc = ac.output.project_acolite_netcdf(l2r_nc)
        proj_time = time.time() - t1

        assert projected_nc is not None, "Reprojection failed"

        # Export to GeoTIFF
        t2 = time.time()
        from osgeo import gdal
        geotiff_dir = out_dir / "geotiff"
        geotiff_dir.mkdir(exist_ok=True)

        ds = gdal.Open(projected_nc)
        subdatasets = ds.GetSubDatasets()
        ds = None

        tif_files = []
        for sds_path, _ in subdatasets:
            varname = sds_path.split(":")[-1]
            if varname in ["x", "y", "transverse_mercator", "lat", "lon"]:
                continue
            outfile = str(geotiff_dir / f"{varname}.tif")
            src_ds = gdal.Open(sds_path)
            if src_ds is None:
                continue
            band = src_ds.GetRasterBand(1)
            nodata = band.GetNoDataValue()
            warp_opts = gdal.WarpOptions(
                format="GTiff",
                srcNodata=nodata,
                dstNodata=float("nan"),
                resampleAlg=gdal.GRA_NearestNeighbour,
                creationOptions=["COMPRESS=DEFLATE", "TILED=YES"],
            )
            gdal.Warp(outfile, src_ds, options=warp_opts)
            src_ds = None
            tif_files.append(outfile)
        export_time = time.time() - t2

        total_time = time.time() - t0
        peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

        return {
            "l2r_nc": l2r_nc,
            "projected_nc": projected_nc,
            "tif_files": tif_files,
            "ac_time": ac_time,
            "proj_time": proj_time,
            "export_time": export_time,
            "total_time": total_time,
            "peak_rss_mb": peak_mb,
        }

    def test_l2r_produced(self, processing_result):
        """ACOLITE produces an L2R NetCDF."""
        assert os.path.exists(processing_result["l2r_nc"])

    def test_reprojection_produced(self, processing_result):
        """Reprojection produces projected NetCDF."""
        assert os.path.exists(processing_result["projected_nc"])

    def test_geotiffs_produced(self, processing_result):
        """Export produces sufficient GeoTIFF files."""
        tifs = processing_result["tif_files"]
        assert len(tifs) >= MIN_COG_FILES, \
            f"Expected >= {MIN_COG_FILES} COGs, got {len(tifs)}"

    def test_geotiffs_not_empty(self, processing_result):
        """Each GeoTIFF is non-trivial in size."""
        for tif in processing_result["tif_files"]:
            size_kb = os.path.getsize(tif) / 1024
            assert size_kb > MIN_COG_SIZE_KB, \
                f"{os.path.basename(tif)} is only {size_kb:.0f} KB"

    def test_geotiff_has_valid_data(self, processing_result):
        """At least one rhos band has non-NaN data."""
        from osgeo import gdal
        rhos_tifs = [t for t in processing_result["tif_files"] if "rhos_" in t]
        assert rhos_tifs, "No rhos GeoTIFF files found"

        ds = gdal.Open(rhos_tifs[0])
        data = ds.GetRasterBand(1).ReadAsArray()
        ds = None

        valid = np.isfinite(data).sum()
        assert valid > 0, "rhos band has no valid (non-NaN) pixels"
        print(f"\n  rhos valid pixels: {valid} / {data.size} ({100*valid/data.size:.1f}%)")

    def test_geotiff_has_projection(self, processing_result):
        """GeoTIFFs have valid CRS (EPSG:32754)."""
        from osgeo import gdal, osr
        rhos_tifs = [t for t in processing_result["tif_files"] if "rhos_" in t]
        ds = gdal.Open(rhos_tifs[0])
        srs = osr.SpatialReference(wkt=ds.GetProjection())
        epsg = srs.GetAuthorityCode(None)
        ds = None
        assert epsg == "32754", f"Expected EPSG:32754, got {epsg}"

    def test_processing_time(self, processing_result):
        """Total processing within time budget."""
        total = processing_result["total_time"]
        print(f"\n  Timing breakdown:")
        print(f"    AC:          {processing_result['ac_time']:.1f}s")
        print(f"    Reproject:   {processing_result['proj_time']:.1f}s")
        print(f"    GeoTIFF:     {processing_result['export_time']:.1f}s")
        print(f"    Total:       {total:.1f}s")
        assert total < MAX_PROCESSING_TIME_S, \
            f"Processing took {total:.0f}s, budget is {MAX_PROCESSING_TIME_S}s"

    def test_memory_usage(self, processing_result):
        """Peak RSS within memory budget."""
        peak = processing_result["peak_rss_mb"]
        print(f"\n  Peak RSS: {peak:.0f} MB (budget: {MAX_PEAK_RSS_MB} MB)")
        assert peak < MAX_PEAK_RSS_MB, \
            f"Peak RSS {peak:.0f} MB exceeds {MAX_PEAK_RSS_MB} MB budget"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--timeout=900", "-s"])
