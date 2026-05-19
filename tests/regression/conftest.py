"""Pytest configuration and shared fixtures for regression testing."""

import pytest
import os
import subprocess
import tempfile
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False,
        help="Run slow tests that require real satellite data"
    )
    parser.addoption(
        "--pace-file", action="store", default=None,
        help="Path to a real PACE OCI L1B NetCDF file for regression"
    )
    parser.addoption(
        "--landsat-file", action="store", default=None,
        help="Path to a Landsat 8/9 bundle directory for regression"
    )
    parser.addoption(
        "--s2-file", action="store", default=None,
        help="Path to a Sentinel-2 .SAFE directory for regression"
    )
    parser.addoption(
        "--tolerance", action="store", default="1e-4", type=float,
        help="Max absolute difference allowed between Python and Rust outputs"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: requires real satellite data (use --runslow)")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        skip = pytest.mark.skip(reason="needs --runslow option")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip)


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def rust_binary(repo_root):
    """Build the Rust binary once per session and return its path."""
    cargo_toml = repo_root / "Cargo.toml"
    if not cargo_toml.exists():
        pytest.skip("No Cargo.toml found — Rust crate not present")

    result = subprocess.run(
        ["cargo", "build", "--release", "--features", "netcdf"],
        cwd=str(repo_root), capture_output=True, text=True, timeout=600,
    )
    if result.returncode != 0:
        pytest.skip(f"Cargo build failed: {result.stderr[:500]}")

    binary = repo_root / "target" / "release" / "acolite-rs"
    if not binary.exists():
        pytest.skip("Rust binary not found after build")
    return binary


@pytest.fixture(scope="session")
def tolerance(request):
    return request.config.getoption("--tolerance")


@pytest.fixture
def tmp_output():
    with tempfile.TemporaryDirectory(prefix="acolite_regtest_") as d:
        yield Path(d)


@pytest.fixture(scope="session")
def pace_file(request):
    path = request.config.getoption("--pace-file")
    if path and Path(path).exists():
        return Path(path)
    # Check env var fallback
    env_path = os.environ.get("ACOLITE_PACE_TEST_FILE")
    if env_path and Path(env_path).exists():
        return Path(env_path)
    return None


@pytest.fixture(scope="session")
def landsat_dir(request):
    path = request.config.getoption("--landsat-file")
    if path and Path(path).is_dir():
        return Path(path)
    env_path = os.environ.get("ACOLITE_LANDSAT_TEST_DIR")
    if env_path and Path(env_path).is_dir():
        return Path(env_path)
    return None


@pytest.fixture(scope="session")
def s2_dir(request):
    path = request.config.getoption("--s2-file")
    if path and Path(path).is_dir():
        return Path(path)
    env_path = os.environ.get("ACOLITE_S2_TEST_DIR")
    if env_path and Path(env_path).is_dir():
        return Path(env_path)
    return None
