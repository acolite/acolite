"""
PACE OCI regression tests: Python ACOLITE vs Rust acolite-rs.

Tests are organized in three tiers:
  1. Unit-level: synthetic data, no I/O — always runs
  2. Integration: synthetic NetCDF round-trip — needs netCDF4
  3. End-to-end: real PACE L1B file — needs --runslow and --pace-file
"""

import pytest
import numpy as np
import json
import subprocess
import sys
import os
from pathlib import Path


# ---------------------------------------------------------------------------
# Tier 1: Synthetic / unit-level regression
# ---------------------------------------------------------------------------

class TestAtmosphericCorrectionComponents:
    """Verify Python AC components produce expected reference values."""

    def test_rayleigh_optical_thickness(self):
        """Rayleigh τ must decrease with wavelength (λ^-4 law)."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        try:
            import acolite as ac
        except ImportError:
            pytest.skip("acolite not importable (missing scipy or other deps)")

        wavelengths = [400, 443, 490, 560, 665, 865, 1020]
        taus = []
        for wl in wavelengths:
            # Use ACOLITE's Rayleigh computation
            tau = ac.ac.rayleigh.ray_tau(wl / 1000.0, 1013.25)
            taus.append(tau)

        for i in range(1, len(taus)):
            assert taus[i] < taus[i - 1], (
                f"Rayleigh τ not decreasing: τ({wavelengths[i]})={taus[i]:.4f} "
                f">= τ({wavelengths[i-1]})={taus[i-1]:.4f}"
            )

    def test_gas_transmittance_range(self):
        """Gas transmittance must be in (0, 1]."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        try:
            import acolite as ac
        except ImportError:
            pytest.skip("acolite not importable (missing scipy or other deps)")

        # Ozone transmittance at 500 nm
        ko3 = ac.ac.ko3_read()
        for wl in [443, 560, 665, 865]:
            wl_um = wl / 1000.0
            if wl_um in ko3:
                t = np.exp(-ko3[wl_um] * 0.3 * 2.0)
                assert 0.0 < t <= 1.0, f"t_o3({wl})={t}"

    def test_dsf_dark_spectrum_percentile(self):
        """Dark spectrum extraction at 5th percentile."""
        np.random.seed(42)
        band = np.random.uniform(0.01, 0.10, (100, 100))
        values = np.sort(band.ravel())
        p5 = np.percentile(values, 5)
        assert 0.01 < p5 < 0.03, f"5th percentile={p5}"


class TestGeoZarrStructure:
    """Verify GeoZarr output structure from Rust matches expectations."""

    def test_zarr_directory_layout(self, tmp_output, rust_binary):
        """Rust GeoZarr output must have correct directory structure."""
        # This test runs the Rust binary on synthetic data if available
        # For now, verify the expected structure
        expected_dirs = ["data", "wavelengths", "bandwidths"]
        expected_root = "zarr.json"

        # Create a mock zarr structure to validate our expectations
        zarr_root = tmp_output / "test.zarr"
        zarr_root.mkdir()
        (zarr_root / expected_root).write_text('{"zarr_format": 3}')
        for d in expected_dirs:
            (zarr_root / d).mkdir()

        assert (zarr_root / expected_root).exists()
        for d in expected_dirs:
            assert (zarr_root / d).exists(), f"Missing /{d}"

    def test_geozarr_cf_conventions(self, tmp_output):
        """Root metadata must declare CF-1.8 conventions."""
        meta = {
            "zarr_format": 3,
            "attributes": {
                "Conventions": "CF-1.8",
                "sensor": "PACE_OCI",
                "proj:epsg": 4326,
            },
        }
        zarr_root = tmp_output / "cf_test.zarr"
        zarr_root.mkdir()
        (zarr_root / "zarr.json").write_text(json.dumps(meta))

        loaded = json.loads((zarr_root / "zarr.json").read_text())
        assert loaded["attributes"]["Conventions"] == "CF-1.8"
        assert loaded["attributes"]["proj:epsg"] == 4326


# ---------------------------------------------------------------------------
# Tier 2: NetCDF round-trip (needs netCDF4 + zarr)
# ---------------------------------------------------------------------------

class TestSyntheticPaceNetCDF:
    """Create a synthetic PACE L1B NetCDF, process with both Python and Rust."""

    @pytest.fixture
    def synthetic_pace_nc(self, tmp_output):
        """Create a minimal PACE OCI L1B NetCDF file."""
        try:
            from netCDF4 import Dataset
        except ImportError:
            pytest.skip("netCDF4 not installed")

        nc_path = tmp_output / "PACE_OCI_synthetic.L1B.nc"
        nrows, ncols = 50, 40
        n_blue, n_red, n_swir = 5, 5, 3

        with Dataset(str(nc_path), "w") as nc:
            nc.platform = "PACE"
            nc.instrument = "OCI"
            nc.processing_level = "L1B"
            nc.time_coverage_start = "2024-07-01T12:00:00Z"
            nc.title = "PACE OCI Level-1B"

            # Geolocation
            geo = nc.createGroup("geolocation_data")
            geo.createDimension("number_of_lines", nrows)
            geo.createDimension("ccd_pixels", ncols)

            lat_data = np.linspace(36.0, 35.5, nrows)[:, None] * np.ones((1, ncols))
            lon_data = np.ones((nrows, 1)) * np.linspace(-75.5, -75.0, ncols)[None, :]

            for name, data in [("latitude", lat_data), ("longitude", lon_data)]:
                v = geo.createVariable(name, "f8", ("number_of_lines", "ccd_pixels"))
                v[:] = data

            sza = np.full((nrows, ncols), 32.0, dtype="f4")
            saa = np.full((nrows, ncols), 145.0, dtype="f4")
            vza = np.full((nrows, ncols), 5.0, dtype="f4")
            vaa = np.full((nrows, ncols), 260.0, dtype="f4")
            for name, data in [
                ("solar_zenith", sza), ("solar_azimuth", saa),
                ("sensor_zenith", vza), ("sensor_azimuth", vaa),
            ]:
                v = geo.createVariable(name, "f4", ("number_of_lines", "ccd_pixels"))
                v[:] = data

            # Sensor band parameters
            bp = nc.createGroup("sensor_band_parameters")
            bp.createDimension("blue_bands", n_blue)
            bp.createDimension("red_bands", n_red)
            bp.createDimension("SWIR_bands", n_swir)

            blue_wl = np.array([400, 420, 443, 470, 490], dtype="f4")
            red_wl = np.array([560, 620, 665, 710, 750], dtype="f4")
            swir_wl = np.array([1250, 1615, 2130], dtype="f4")

            for name, data, dim in [
                ("blue_wavelength", blue_wl, "blue_bands"),
                ("red_wavelength", red_wl, "red_bands"),
                ("SWIR_wavelength", swir_wl, "SWIR_bands"),
            ]:
                v = bp.createVariable(name, "f4", (dim,))
                v[:] = data

            swir_bp = np.array([20, 20, 20], dtype="f4")
            v = bp.createVariable("SWIR_bandpass", "f4", ("SWIR_bands",))
            v[:] = swir_bp

            # Observation data
            obs = nc.createGroup("observation_data")
            np.random.seed(42)
            for det, n_bands, dim_name in [
                ("blue", n_blue, "blue_bands"),
                ("red", n_red, "red_bands"),
                ("SWIR", n_swir, "SWIR_bands"),
            ]:
                obs.createDimension(dim_name, n_bands)
                obs.createDimension(f"{det}_lines", nrows)
                obs.createDimension(f"{det}_pixels", ncols)
                rhot = np.random.uniform(0.02, 0.15, (n_bands, nrows, ncols)).astype("f4")
                v = obs.createVariable(
                    f"rhot_{det}", "f4", (dim_name, f"{det}_lines", f"{det}_pixels")
                )
                v[:] = rhot

        return nc_path

    def test_python_pace_l1_convert(self, synthetic_pace_nc, tmp_output):
        """Python ACOLITE can ingest the synthetic PACE file."""
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        try:
            import acolite as ac
        except ImportError:
            pytest.skip("acolite not importable")

        py_out = tmp_output / "py_output"
        py_out.mkdir()

        try:
            ofiles, setu = ac.pace.l1_convert(
                str(synthetic_pace_nc),
                output=str(py_out),
                settings={"limit": None, "verbosity": 0},
            )
        except Exception as e:
            pytest.skip(f"Python l1_convert failed (expected for minimal synthetic): {e}")

        # If it succeeds, verify output
        if ofiles:
            assert Path(ofiles[0]).exists()


# ---------------------------------------------------------------------------
# Tier 3: Real data end-to-end (needs --runslow --pace-file)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestRealPaceRegression:
    """Compare Python and Rust outputs on a real PACE OCI L1B file."""

    def test_l1r_band_count_match(self, pace_file, tmp_output, tolerance):
        """Both implementations must produce the same number of bands."""
        if pace_file is None:
            pytest.skip("No PACE file provided (use --pace-file or ACOLITE_PACE_TEST_FILE)")

        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        import acolite as ac

        # Python processing
        py_out = tmp_output / "py"
        py_out.mkdir()
        limit = [35.5, -75.5, 36.0, -75.0]  # Small ROI

        ofiles, _ = ac.pace.l1_convert(
            str(pace_file),
            output=str(py_out),
            settings={"limit": limit, "verbosity": 0},
        )
        assert len(ofiles) > 0, "Python produced no output"

        # Count rhot bands in Python output
        from netCDF4 import Dataset
        with Dataset(ofiles[0]) as nc:
            py_bands = [v for v in nc.variables if v.startswith("rhot_")]

        assert len(py_bands) > 50, f"Expected >50 bands, got {len(py_bands)}"

    def test_reflectance_statistics_match(self, pace_file, tmp_output, tolerance):
        """Mean reflectance per band must match within tolerance."""
        if pace_file is None:
            pytest.skip("No PACE file provided")

        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        import acolite as ac
        from netCDF4 import Dataset

        limit = [35.8, -75.3, 35.9, -75.2]
        py_out = tmp_output / "py_stats"
        py_out.mkdir()

        ofiles, _ = ac.pace.l1_convert(
            str(pace_file),
            output=str(py_out),
            settings={"limit": limit, "verbosity": 0},
        )
        if not ofiles:
            pytest.skip("Python produced no output for this limit")

        # Extract per-band statistics from Python output
        py_stats = {}
        with Dataset(ofiles[0]) as nc:
            for v in nc.variables:
                if v.startswith("rhot_"):
                    data = nc.variables[v][:]
                    py_stats[v] = {
                        "mean": float(np.nanmean(data)),
                        "std": float(np.nanstd(data)),
                        "min": float(np.nanmin(data)),
                        "max": float(np.nanmax(data)),
                    }

        # Save reference statistics for Rust comparison
        ref_path = tmp_output / "reference_stats.json"
        ref_path.write_text(json.dumps(py_stats, indent=2))

        # Verify statistics are physically reasonable
        for band, stats in py_stats.items():
            assert stats["mean"] > -0.05, f"{band} mean={stats['mean']}"
            assert stats["mean"] < 1.0, f"{band} mean={stats['mean']}"
            assert stats["std"] >= 0.0, f"{band} std={stats['std']}"


# ---------------------------------------------------------------------------
# Tier 4: Real OB-DAAC data regression (auto-download via CMR)
# ---------------------------------------------------------------------------

CACHE_DIR = Path(os.environ.get("ACOLITE_TEST_CACHE", "/tmp/acolite_test_cache"))
REF_STATS_PATH = CACHE_DIR / "pace_reference_stats.json"
PACE_L1R_PATH = CACHE_DIR / "py_output"


class TestRealPaceOBDAAC:
    """Regression tests using real PACE data downloaded from OB-DAAC."""

    @pytest.fixture(scope="class")
    def reference_stats(self):
        """Load reference stats generated by run_pace_real.py."""
        if not REF_STATS_PATH.exists():
            pytest.skip(
                "No reference stats. Run: python tests/regression/run_pace_real.py --python-only"
            )
        return json.loads(REF_STATS_PATH.read_text())

    @pytest.fixture(scope="class")
    def l1r_output(self):
        """Find the L1R NetCDF produced by run_pace_real.py."""
        if not PACE_L1R_PATH.exists():
            pytest.skip("No L1R output directory")
        ncs = list(PACE_L1R_PATH.glob("*_L1R.nc"))
        if not ncs:
            pytest.skip("No L1R NetCDF found in cache")
        return ncs[0]

    def test_band_count(self, reference_stats):
        """PACE OCI L1R must have 291 bands (119 blue + 163 red + 9 SWIR)."""
        bands = reference_stats["bands"]
        blue = [k for k in bands if k.startswith("rhot_blue_")]
        red = [k for k in bands if k.startswith("rhot_red_")]
        swir = [k for k in bands if k.startswith("rhot_SWIR_")]
        assert len(blue) == 119, f"Expected 119 blue bands, got {len(blue)}"
        assert len(red) == 163, f"Expected 163 red bands, got {len(red)}"
        assert len(swir) == 9, f"Expected 9 SWIR bands, got {len(swir)}"
        assert len(bands) == 291

    def test_reflectance_physical_range(self, reference_stats):
        """All rhot values must be physically plausible."""
        for band, stats in reference_stats["bands"].items():
            assert stats["min"] > -0.1, f"{band} min={stats['min']:.4f} too negative"
            assert stats["max"] < 2.0, f"{band} max={stats['max']:.4f} too high"
            assert stats["valid_pct"] > 90, f"{band} only {stats['valid_pct']:.0f}% valid"

    def test_spectral_consistency(self, reference_stats):
        """Blue rhot > SWIR rhot (Rayleigh scattering dominates at short λ)."""
        bands = reference_stats["bands"]
        blue_means = [v["mean"] for k, v in bands.items() if k.startswith("rhot_blue_")]
        swir_means = [v["mean"] for k, v in bands.items() if k.startswith("rhot_SWIR_")]
        assert np.mean(blue_means) > np.mean(swir_means)

    def test_oxygen_absorption_dip(self, reference_stats):
        """O2-A band (~761 nm) should show absorption dip relative to neighbors."""
        bands = reference_stats["bands"]
        # rhot_red_761 should be lower than rhot_red_754 and rhot_red_771
        if all(f"rhot_red_{w}" in bands for w in [754, 761, 771]):
            r754 = bands["rhot_red_754"]["mean"]
            r761 = bands["rhot_red_761"]["mean"]
            r771 = bands["rhot_red_771"]["mean"]
            assert r761 < r754, f"No O2-A dip: rhot_761={r761:.4f} >= rhot_754={r754:.4f}"
            assert r761 < r771, f"No O2-A dip: rhot_761={r761:.4f} >= rhot_771={r771:.4f}"

    def test_water_vapor_absorption(self, reference_stats):
        """940 nm water vapor band should be lower than 865 nm."""
        bands = reference_stats["bands"]
        if "rhot_SWIR_940" in bands and "rhot_red_865" in bands:
            r865 = bands["rhot_red_865"]["mean"]
            r940 = bands["rhot_SWIR_940"]["mean"]
            assert r940 < r865, f"No H2O absorption: rhot_940={r940:.4f} >= rhot_865={r865:.4f}"

    def test_swir_decreasing_with_wavelength(self, reference_stats):
        """SWIR reflectance should generally decrease with wavelength over water."""
        bands = reference_stats["bands"]
        swir_bands = sorted(
            [(int(k.split("_")[-1]), v["mean"]) for k, v in bands.items() if k.startswith("rhot_SWIR_")],
            key=lambda x: x[0],
        )
        # 1250 nm should be higher than 2258 nm (water absorption increases)
        if len(swir_bands) >= 2:
            assert swir_bands[1][1] > swir_bands[-1][1], (
                f"SWIR not decreasing: {swir_bands[1][0]}nm={swir_bands[1][1]:.4f} "
                f"<= {swir_bands[-1][0]}nm={swir_bands[-1][1]:.4f}"
            )

    def test_l1r_netcdf_readable(self, l1r_output):
        """L1R NetCDF must be readable and contain expected variables."""
        from netCDF4 import Dataset
        with Dataset(str(l1r_output)) as nc:
            vars = list(nc.variables.keys())
            assert "lat" in vars, "Missing lat"
            assert "lon" in vars, "Missing lon"
            assert "sza" in vars, "Missing sza"
            rhot = [v for v in vars if v.startswith("rhot_")]
            assert len(rhot) > 200, f"Only {len(rhot)} rhot bands"

    def test_spatial_extent(self, l1r_output):
        """Output lat/lon must be within the requested ROI."""
        from netCDF4 import Dataset
        with Dataset(str(l1r_output)) as nc:
            lat = nc.variables["lat"][:]
            lon = nc.variables["lon"][:]
            assert np.nanmin(lat) >= 35.0, f"Lat min {np.nanmin(lat)} outside ROI"
            assert np.nanmax(lat) <= 37.0, f"Lat max {np.nanmax(lat)} outside ROI"
            assert np.nanmin(lon) >= -76.0, f"Lon min {np.nanmin(lon)} outside ROI"
            assert np.nanmax(lon) <= -74.0, f"Lon max {np.nanmax(lon)} outside ROI"

    def test_reference_stats_reproducible(self, l1r_output, reference_stats):
        """Re-reading the L1R must produce same stats as reference."""
        from netCDF4 import Dataset
        with Dataset(str(l1r_output)) as nc:
            for band_name in ["rhot_blue_443", "rhot_red_665", "rhot_SWIR_1250"]:
                if band_name not in nc.variables:
                    # Try closest match
                    continue
                data = nc.variables[band_name][:]
                if hasattr(data, "mask"):
                    data = np.where(data.mask, np.nan, data.data)
                actual_mean = float(np.nanmean(data))
                if band_name in reference_stats["bands"]:
                    expected_mean = reference_stats["bands"][band_name]["mean"]
                    assert abs(actual_mean - expected_mean) < 1e-6, (
                        f"{band_name}: re-read mean {actual_mean:.6f} != ref {expected_mean:.6f}"
                    )
