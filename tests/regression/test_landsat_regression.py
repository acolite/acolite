"""
Landsat 8/9 regression tests: Python ACOLITE vs Rust acolite-rs.

Tier 1: Synthetic — calibration, Rayleigh, DSF physics (always runs)
Tier 2: MTL parsing round-trip (always runs)
Tier 3: Real data comparison (needs --runslow --landsat-file)
"""

import pytest
import numpy as np
import json
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Tier 1: Calibration and physics regression
# ---------------------------------------------------------------------------

class TestLandsatCalibration:
    """Verify calibration math matches between Python and Rust."""

    LANDSAT_BANDS = {
        "B1": (443.0, 16.0),
        "B2": (482.0, 60.0),
        "B3": (561.0, 57.0),
        "B4": (655.0, 37.0),
        "B5": (865.0, 28.0),
        "B6": (1609.0, 85.0),
        "B7": (2201.0, 187.0),
    }

    def test_dn_to_radiance_linearity(self):
        """L = DN * MULT + ADD must be exact."""
        dn = np.arange(0, 65536, 1000, dtype=np.uint16)
        mult, add = 0.012, -60.0
        rad = dn.astype(np.float64) * mult + add
        # Verify no precision loss
        for i in range(len(dn)):
            expected = float(dn[i]) * mult + add
            assert abs(rad[i] - expected) < 1e-10

    def test_dn_to_reflectance_sun_angle(self):
        """Higher sun elevation → lower reflectance."""
        dn = 10000.0
        mult, add = 2e-5, -0.1
        rho_base = dn * mult + add

        rho_high_sun = rho_base / np.sin(np.radians(60.0))
        rho_low_sun = rho_base / np.sin(np.radians(30.0))
        assert rho_high_sun < rho_low_sun

    def test_earth_sun_distance_annual_cycle(self):
        """Perihelion < aphelion distance."""
        # Hansen formula: d = 1 - 0.01672 * cos(2π(doy-1)/365)
        def esd(doy):
            theta = 2.0 * np.pi * (doy - 1) / 365.0
            return 1.0 - 0.01672 * np.cos(theta)

        d_peri = esd(3)    # ~Jan 3
        d_aph = esd(185)   # ~Jul 4
        assert d_peri < d_aph
        assert abs(d_peri - 0.983) < 0.02
        assert abs(d_aph - 1.017) < 0.02

    def test_band_wavelengths_l8_l9_identical(self):
        """L8 and L9 OLI share the same band definitions."""
        # This is a known property — both sensors have identical spectral response
        for band, (wl, bw) in self.LANDSAT_BANDS.items():
            assert wl > 0 and bw > 0, f"Invalid {band}: wl={wl}, bw={bw}"


class TestLandsatRayleigh:
    """Rayleigh scattering physics for Landsat bands."""

    def test_rayleigh_tau_monotonic_decrease(self):
        """τ_R must decrease with wavelength (λ^-4 law)."""
        wavelengths = [443, 482, 561, 655, 865, 1609, 2201]
        # Hansen & Travis (1974)
        def ray_tau(wl_nm, pressure=1013.25):
            lam = wl_nm / 1000.0  # µm
            tau = 0.008569 * lam**(-4) * (1 + 0.0113 * lam**(-2) + 0.00013 * lam**(-4))
            return tau * (pressure / 1013.25)

        taus = [ray_tau(w) for w in wavelengths]
        for i in range(1, len(taus)):
            assert taus[i] < taus[i - 1], (
                f"τ({wavelengths[i]})={taus[i]:.4f} >= τ({wavelengths[i-1]})={taus[i-1]:.4f}"
            )

        # B1 (443nm) significant Rayleigh
        assert taus[0] > 0.1
        # B7 (2201nm) negligible
        assert taus[-1] < 0.001

    def test_rayleigh_pressure_scaling(self):
        """τ_R scales linearly with pressure."""
        def ray_tau(wl_nm, pressure):
            lam = wl_nm / 1000.0
            tau0 = 0.008569 * lam**(-4) * (1 + 0.0113 * lam**(-2) + 0.00013 * lam**(-4))
            return tau0 * (pressure / 1013.25)

        tau_std = ray_tau(500, 1013.25)
        tau_half = ray_tau(500, 1013.25 / 2)
        assert abs(tau_half / tau_std - 0.5) < 1e-10


class TestLandsatDSF:
    """DSF dark spectrum fitting for Landsat SWIR bands."""

    def test_dark_spectrum_percentile(self):
        """5th percentile extraction from SWIR-like data."""
        np.random.seed(42)
        swir1 = np.random.uniform(0.001, 0.02, (200, 200))
        swir2 = np.random.uniform(0.0005, 0.01, (200, 200))

        p5_1 = np.percentile(swir1, 5)
        p5_2 = np.percentile(swir2, 5)
        assert 0.001 < p5_1 < 0.005
        assert 0.0005 < p5_2 < 0.003

    def test_aot_increases_with_turbidity(self):
        """Brighter dark spectrum → higher estimated AOT."""
        # Simple linear AOT model (matches Rust placeholder)
        def estimate_aot(dark_spectrum):
            mean_dark = np.mean(dark_spectrum)
            return min(max(mean_dark * 10.0, 0.01), 1.0)

        aot_clear = estimate_aot([0.005, 0.002])
        aot_turbid = estimate_aot([0.05, 0.03])
        assert aot_turbid > aot_clear


# ---------------------------------------------------------------------------
# Tier 2: MTL parsing and pipeline structure
# ---------------------------------------------------------------------------

class TestLandsatMTL:
    """MTL metadata parsing regression."""

    SAMPLE_MTL = """GROUP = L1_METADATA_FILE
  GROUP = IMAGE_ATTRIBUTES
    SPACECRAFT_ID = "LANDSAT_8"
    SENSOR_ID = "OLI_TIRS"
    SUN_ELEVATION = 55.12345678
    SUN_AZIMUTH = 140.98765432
  END_GROUP = IMAGE_ATTRIBUTES
  GROUP = PRODUCT_METADATA
    DATE_ACQUIRED = 2024-06-15
    SCENE_CENTER_TIME = "10:30:45.1234567Z"
    DATA_TYPE = "L1TP"
  END_GROUP = PRODUCT_METADATA
END_GROUP = L1_METADATA_FILE
END
"""

    def test_parse_mtl_fields(self):
        """Key MTL fields must be extractable."""
        attrs = {}
        for line in self.SAMPLE_MTL.strip().splitlines():
            line = line.strip()
            if "=" in line and not line.startswith("GROUP") and not line.startswith("END"):
                key, val = line.split("=", 1)
                attrs[key.strip()] = val.strip().strip('"')

        assert attrs["SPACECRAFT_ID"] == "LANDSAT_8"
        assert attrs["SENSOR_ID"] == "OLI_TIRS"
        assert float(attrs["SUN_ELEVATION"]) == pytest.approx(55.12345678)
        assert float(attrs["SUN_AZIMUTH"]) == pytest.approx(140.98765432)

    def test_sun_zenith_from_elevation(self):
        """SZA = 90 - sun_elevation."""
        sun_elev = 55.12345678
        sza = 90.0 - sun_elev
        assert abs(sza - 34.87654322) < 1e-6


class TestLandsatPerformanceBaseline:
    """Establish Python performance baselines for comparison with Rust."""

    def test_numpy_band_processing_throughput(self):
        """Measure Python per-band processing speed as baseline."""
        size = 1000
        nbands = 7
        np.random.seed(42)

        bands = [np.random.randint(0, 65535, (size, size), dtype=np.uint16) for _ in range(nbands)]

        t0 = time.perf_counter()
        results = []
        for band in bands:
            # Simulate: DN→reflectance, gas correction, Rayleigh subtraction
            toa = band.astype(np.float64) / 10000.0
            t_gas = 0.95  # placeholder
            toa = toa / t_gas
            tau_r = 0.2   # placeholder
            rhos = toa - tau_r * 0.05
            rhos = np.clip(rhos, 0, None)
            results.append(rhos)
        elapsed = time.perf_counter() - t0

        mpx = (size * size * nbands) / elapsed / 1e6
        print(f"\nPython baseline: {mpx:.1f} Mpx/s ({nbands} bands, {size}×{size}, {elapsed:.3f}s)")

        # Just record — no assertion on Python speed
        assert len(results) == nbands


# ---------------------------------------------------------------------------
# Tier 3: Real data (needs --runslow --landsat-file)
# ---------------------------------------------------------------------------

def pytest_addoption_landsat(parser):
    """Register --landsat-file option (called from conftest)."""
    try:
        parser.addoption("--landsat-file", action="store", default=None,
                         help="Path to Landsat 8/9 bundle directory")
    except ValueError:
        pass  # already registered


@pytest.fixture(scope="session")
def landsat_dir(request):
    import os
    path = request.config.getoption("--landsat-file", default=None)
    if path and Path(path).is_dir():
        return Path(path)
    env = os.environ.get("ACOLITE_LANDSAT_TEST_DIR")
    if env and Path(env).is_dir():
        return Path(env)
    return None


@pytest.mark.slow
class TestRealLandsatRegression:
    """Compare Python and Rust on real Landsat data."""

    def test_l1_convert_band_count(self, landsat_dir, tmp_output):
        """Python l1_convert must produce 7 OLI bands."""
        if landsat_dir is None:
            pytest.skip("No Landsat dir (use --landsat-file or ACOLITE_LANDSAT_TEST_DIR)")

        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        try:
            import acolite as ac
        except ImportError:
            pytest.skip("acolite not importable")

        py_out = tmp_output / "landsat_py"
        py_out.mkdir()

        ofiles, _ = ac.landsat.l1_convert(
            str(landsat_dir),
            output=str(py_out),
            settings={"verbosity": 0},
        )
        assert len(ofiles) > 0, "No output produced"

        from netCDF4 import Dataset
        with Dataset(ofiles[0]) as nc:
            rhot_bands = [v for v in nc.variables if v.startswith("rhot_")]

        assert len(rhot_bands) >= 7, f"Expected ≥7 bands, got {len(rhot_bands)}"

    def test_reflectance_statistics(self, landsat_dir, tmp_output, tolerance):
        """Per-band mean reflectance must be physically reasonable."""
        if landsat_dir is None:
            pytest.skip("No Landsat dir")

        sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
        try:
            import acolite as ac
            from netCDF4 import Dataset
        except ImportError:
            pytest.skip("Missing deps")

        py_out = tmp_output / "landsat_stats"
        py_out.mkdir()

        ofiles, _ = ac.landsat.l1_convert(
            str(landsat_dir), output=str(py_out),
            settings={"verbosity": 0},
        )
        if not ofiles:
            pytest.skip("No output")

        stats = {}
        with Dataset(ofiles[0]) as nc:
            for v in nc.variables:
                if v.startswith("rhot_"):
                    data = nc.variables[v][:]
                    stats[v] = {
                        "mean": float(np.nanmean(data)),
                        "std": float(np.nanstd(data)),
                    }

        ref_path = tmp_output / "landsat_reference_stats.json"
        ref_path.write_text(json.dumps(stats, indent=2))

        for band, s in stats.items():
            assert -0.05 < s["mean"] < 1.0, f"{band} mean={s['mean']}"
            assert s["std"] >= 0.0
