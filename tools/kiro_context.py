#!/usr/bin/env python3
"""
kiro_context.py — Compact Kiro resume prompt generator.

Introspects the repo to auto-detect implementation state and emits a
dense context string that fits inside a single ACP session/prompt without
hitting context-window limits.  Use this whenever a Kiro session is
starting fresh or has lost history.

Usage
-----
    # Print compact context to stdout
    python tools/kiro_context.py

    # Target the next priority task explicitly
    python tools/kiro_context.py --next s3-aot-fix
    python tools/kiro_context.py --next l2r-output
    python tools/kiro_context.py --next ancillary-interp
    python tools/kiro_context.py --next dem-pressure
    python tools/kiro_context.py --next settings-file
    python tools/kiro_context.py --next l2w
    python tools/kiro_context.py --next emit
    python tools/kiro_context.py --next perf-ci

    # Pipe directly into the agent harness
    python tools/agent_harness.py \\
        --task "$(python tools/kiro_context.py --next sentinel3)" \\
        --kiro-agent rust-developer

    # Auto-resume: let the script pick the highest-priority incomplete task
    python tools/agent_harness.py \\
        --task "$(python tools/kiro_context.py --auto)" \\
        --kiro-agent rust-developer \\
        --auto-approve
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# State detection — derived from actual files, not from the roadmap docs
# ---------------------------------------------------------------------------

def _file_exists(*parts: str) -> bool:
    return (ROOT / Path(*parts)).exists()


def _line_count(path: str) -> int:
    p = ROOT / path
    if not p.exists():
        return 0
    return sum(1 for _ in open(p))


def _grep_count(pattern: str, *paths: str) -> int:
    """Count lines matching pattern across one or more file globs."""
    total = 0
    for glob in paths:
        for f in ROOT.glob(glob):
            try:
                total += sum(1 for l in open(f) if re.search(pattern, l))
            except Exception:
                pass
    return total


def detect_state() -> dict:
    """Return a dict describing what is done and what is pending."""

    loaders = {
        "landsat":   _file_exists("src/loader/landsat.rs"),
        "sentinel2": _file_exists("src/loader/sentinel2.rs"),
        "pace":      _file_exists("src/loader/pace.rs"),
        "sentinel3": _file_exists("src/loader/sentinel3.rs"),
        "emit":      _file_exists("src/loader/emit.rs"),
        "prisma":    _file_exists("src/loader/prisma.rs"),
    }

    # sentinel3 loader check: stub sensor def only when file is tiny
    s3_loader_real = loaders["sentinel3"] and _line_count("src/loader/sentinel3.rs") > 30
    loaders["sentinel3"] = s3_loader_real

    ac_modules = {
        "interp":       _file_exists("src/ac/interp.rs"),
        "aerlut":       _file_exists("src/ac/aerlut.rs"),
        "gas_lut":      _file_exists("src/ac/gas_lut.rs"),
        "dsf":          _file_exists("src/ac/dsf.rs"),
        "rayleigh":     _file_exists("src/ac/rayleigh.rs"),
        "calibration":  _file_exists("src/ac/calibration.rs"),
        "ancillary":    _file_exists("src/ac/ancillary.rs"),
        "glint":        _file_exists("src/ac/glint.rs"),
    }

    ancillary_wired = ac_modules["ancillary"] and bool(
        _grep_count(r"ancillary", "src/pipeline.rs")
    )
    # ancillary spatial interpolation: check for bilinear interp logic in ancillary.rs
    ancillary_interp_done = ac_modules["ancillary"] and bool(
        _grep_count(r"interp_gmao|bilinear|spatial_interp|RegularGrid", "src/ac/ancillary.rs")
    )

    writers = {
        "cog":      _file_exists("src/writer/cog.rs"),
        "geozarr":  _file_exists("src/writer/geozarr.rs"),
        "netcdf_l2r": _file_exists("src/writer/netcdf_l2r.rs"),
    }

    rust_tests   = _grep_count(r"#\[test\]", "src/**/*.rs", "tests/*.rs")
    python_tests = _grep_count(r"def test_", "tests/regression/**/*.py")

    # ROI: all 4 loaders now accept a limit parameter — check each loader file
    roi_wired = all(
        bool(_grep_count(r"\blimit\b", f"src/loader/{sensor}.rs"))
        for sensor in ("landsat", "sentinel2", "pace", "sentinel3")
        if _file_exists(f"src/loader/{sensor}.rs")
    )

    l2w_exists = _file_exists("src/ac/l2w.rs")
    perf_ci    = _file_exists(".github/workflows/perf.yml")
    dem_pressure = bool(_grep_count(r"pressure_from_elevation|pressure_elevation", "src/ac/ancillary.rs"))
    settings_file = bool(_grep_count(r"settings.txt|settings_file|from_settings", "src/main.rs"))
    # S3 AOT gap: check if the RGI fix has landed (compare_rgi test or aot_gap comment)
    s3_aot_fixed = bool(_grep_count(r"aot_gap|compare_rgi|s3_aot", "src/ac/interp.rs", "tests/sentinel3_e2e.rs"))

    return {
        "loaders":               loaders,
        "ac_modules":            ac_modules,
        "ancillary_wired":       ancillary_wired,
        "ancillary_interp_done": ancillary_interp_done,
        "writers":               writers,
        "roi_wired":             roi_wired,
        "l2w_exists":            l2w_exists,
        "perf_ci":               perf_ci,
        "dem_pressure":          dem_pressure,
        "settings_file":         settings_file,
        "s3_aot_fixed":          s3_aot_fixed,
        "rust_tests":            rust_tests,
        "python_tests":          python_tests,
    }


# ---------------------------------------------------------------------------
# Next-task definitions — ordered by roadmap priority
# ---------------------------------------------------------------------------

NEXT_TASKS = {
    "s3-aot-fix": {
        "short": "Fix S3 OLCI AOT gap — push R from 0.983 to >0.999",
        "detail": (
            "S3 OLCI achieves R=0.983, below the R>0.999 target. Root cause: AOT 0.076 (Rust) "
            "vs 0.057 (Python) from differences between src/ac/interp.rs RegularGridInterpolator "
            "and scipy's RGI on the 5-D aerosol LUT. "
            "Investigate: (1) add a pytest test feeding identical grid points to both Rust RGI and "
            "scipy.interpolate.RegularGridInterpolator to isolate the gap; "
            "(2) verify axes are sorted ascending in src/ac/interp.rs (scipy requires this); "
            "(3) check DSF dark-spectrum percentile uses the same pixel count across all bands jointly "
            "(n=200 darkest pixels joint, not per-band); "
            "(4) compare per-band romix and dark spectrum values via a --debug flag in dsf.rs. "
            "Target: R>0.999, RMSE<0.002 for S3 OLCI matching other sensors. "
            "Run: cargo test --test sentinel3_e2e && pytest tests/regression/test_s3_rust_vs_python.py -v -s"
        ),
    },
    "l2r-output": {
        "short": "Add NetCDF L2R output writer matching Python ACOLITE format",
        "detail": (
            "Create src/writer/netcdf_l2r.rs behind #[cfg(feature='full-io')]. "
            "Schema must match Python ACOLITE gem.gem() reader exactly: "
            "global attrs (sensor, isodate, acolite_version, limit), "
            "dimensions (y, x), variables (lat f64, lon f64, rhos_NNN f32, rhot_NNN f32 per band). "
            "Variable naming: rhos_443, rhot_443 etc (integer wavelength in nm). "
            "Compression: zlib level 4, chunked 256x256. "
            "Add write_netcdf_l2r() to src/writer/mod.rs dispatch. "
            "Add --output-format netcdf CLI arg. "
            "Round-trip test: Python netCDF4 opens Rust output and reads rhos_ variables within 1e-4."
        ),
    },
    "ancillary-interp": {
        "short": "Port ancillary spatial+temporal interpolation (interp_gmao.py)",
        "detail": (
            "src/ac/ancillary.rs downloads files but uses scene-level defaults, not spatially "
            "interpolated values. Python uses acolite/ac/ancillary/interp_gmao.py (63 lines): "
            "bilinear RGI on lon/lat grid, linear temporal between two 6-hourly MERRA2 files. "
            "Port to Rust: fn interp_gmao(files: &[PathBuf], lon: f64, lat: f64, isodate: &str) -> Result<Ancillary>. "
            "Steps: (1) read lon/lat grid + time_coverage_start from each GMAO NetCDF; "
            "(2) bilinear spatial interp at (lon, lat) using src/ac/interp.rs RGI on 2D lon/lat grid; "
            "(3) linear temporal interp between two bracketing 6-hourly files. "
            "Variables: PS (pressure hPa), TO3 (ozone DU), TQV (water vapour g/cm2), "
            "U10M + V10M (wind components, compute speed = sqrt(u^2+v^2)). "
            "Wire into ancillary::fetch() in pipeline.rs. "
            "Test: known MERRA2 file -> assert interpolated values within 1% of Python reference."
        ),
    },
    "dem-pressure": {
        "short": "Port DEM-derived surface pressure (pressure_elevation.py)",
        "detail": (
            "Port acolite/ac/pressure_elevation.py (33 lines) to src/ac/ancillary.rs. "
            "Standard atmosphere barometric formula: "
            "p = 101325 * exp(-h / h0) / 100 [hPa] where h0 = (R*T0)/(M*g) = 8434.5 m. "
            "Constants: p0=101325 Pa, T0=288.15 K, g=9.80665 m/s2, M=0.0289644 kg/mol, R=8.31447 J/mol/K. "
            "Add fn pressure_from_elevation(elevation_m: f64) -> f64. "
            "Source elevation from SRTM DEM if available at scene centre lat/lon, else fall back to "
            "GMAO PS field, else 1013.25 hPa default. "
            "Add --dem-pressure CLI flag. "
            "Tests: 0m -> 1013.25 hPa, 1000m -> ~899 hPa, 2000m -> ~795 hPa (4 sig figs)."
        ),
    },
    "settings-file": {
        "short": "Add Python ACOLITE-compatible settings file CLI",
        "detail": (
            "Python ACOLITE is driven by flat key=value settings files. Rust only takes CLI flags. "
            "Add --settings /path/to/settings.txt to src/main.rs. "
            "Parse format: one 'key = value' or 'key=value' per line, # comments. "
            "Map Python keys to Rust ProcessingConfig: inputfile, output, limit (S,W,N,E), "
            "l2w_parameters (list), dsf_aot_compute (min|fixed), glint_correction (bool), "
            "ancillary_data (bool), dem_pressure (bool). "
            "Precedence: CLI flags override settings file. "
            "Test: write minimal settings.txt, run Rust with --settings, assert ProcessingConfig matches. "
            "Add examples/from_settings_file.rs demonstrating the workflow."
        ),
    },
    "l2w": {
        "short": "Port chlorophyll-a OC3 water product (first L2W product)",
        "detail": (
            "Port OC3 chlorophyll-a from acolite/acolite/acolite_l2w.py to src/ac/l2w.rs. "
            "fn chl_oc3(rhos: &HashMap<u32, Array2<f32>>) -> Result<Array2<f32>>. "
            "Band ratio: R = log10(max(rhos[443], rhos[490]) / rhos[560]). "
            "Polynomial: chl = 10^(0.2515 - 2.3798*R + 1.5317*R^2 - 0.1428*R^3 - 0.9822*R^4). "
            "Edge cases: rhos <= 0 or missing bands -> f32::NAN (not an error). "
            "For sensors without exact 443/490/560nm bands, use nearest within 15nm. "
            "Wire via --l2w-parameters chl_oc3 CLI flag. Output as extra band with units='mg m-3'. "
            "Tests: known rhos -> expected Chl within 1e-4; all-NaN input -> all-NaN output (no panic)."
        ),
    },
    "emit": {
        "short": "Port EMIT L1B hyperspectral loader",
        "detail": (
            "Port EMIT to src/loader/emit.rs following src/loader/pace.rs. "
            "File: EMIT_L1B_RAD_*_RDN_*.nc. "
            "Variable: 'radiance' shape (lines, samples, bands) — 285 bands 380-2500nm. "
            "Wavelengths: sensor_band_parameters/wavelengths (1D, 285 values). "
            "FWHM: sensor_band_parameters/fwhm (for RSR Gaussian convolution). "
            "Geometry: location/lat (lines, samples), location/lon, location/elev. "
            "Bad bands: sensor_band_parameters/good_wavelengths (bool array). "
            "DN to TOA reflectance: rhot = (radiance * pi * d^2) / (irradiance * cos(sza)). "
            "Sensor definition: src/sensors/emit.rs, 285 bands, emit_NNN naming. "
            "Download via LP DAAC CMR search (EMITL1BRAD_001) using existing cmr.rs pattern. "
            "Output: GeoZarr V3 (285 bands > 50 threshold). "
            "Add tests/emit_e2e.rs with synthetic fixture."
        ),
    },
    "perf-ci": {
        "short": "Add GitHub Actions performance regression gate",
        "detail": (
            "Create .github/workflows/perf.yml: trigger on push to main and PR. "
            "Steps: checkout -> cargo bench --bench performance 2>&1 | tee bench.txt. "
            "Store baselines in benches/performance_baseline.json: "
            "{\"landsat_mpx_s\": 940000, \"s2_mpx_s\": 576000, \"pace_mpx_s\": 18700000}. "
            "Parse bench.txt, compare to baseline, fail if any metric drops >10%. "
            "On main-branch merge: update baseline JSON and commit via github-actions bot. "
            "Also add memory check: /usr/bin/time -v cargo run 2>&1 | grep 'Maximum resident'."
        ),
    },
}

# Priority order for --auto (highest priority first)
PRIORITY = [
    "s3-aot-fix",       # accuracy gap: R=0.983 below R>0.999 target
    "l2r-output",       # unblocks Python L2W consumption of Rust output
    "ancillary-interp", # spatial+temporal MERRA2 interpolation
    "dem-pressure",     # barometric pressure from SRTM elevation
    "settings-file",    # Python-compatible settings file CLI
    "l2w",              # chlorophyll-a OC3 (first water product)
    "emit",             # EMIT hyperspectral loader
    "perf-ci",          # CI performance regression gate
]


def is_done(key: str, state: dict) -> bool:
    if key == "s3-aot-fix":
        return state["s3_aot_fixed"]
    if key == "l2r-output":
        return state["writers"].get("netcdf_l2r", False)
    if key == "ancillary-interp":
        return state["ancillary_interp_done"]
    if key == "dem-pressure":
        return state["dem_pressure"]
    if key == "settings-file":
        return state["settings_file"]
    if key == "l2w":
        return state["l2w_exists"]
    if key == "emit":
        return state["loaders"].get("emit", False)
    if key == "perf-ci":
        return state["perf_ci"]
    return False


# ---------------------------------------------------------------------------
# Compact context block
# ---------------------------------------------------------------------------

PHYSICS = (
    "rhos=(rhot/tt_gas-romix)/(dutott+astot*(rhot/tt_gas-romix)) | "
    "DSF=intercept(200px darkest) | wave_range=(400,900)nm | tgas_min=0.85 | "
    "AOT=min-RMSD model selection | NaN when ds<romix(tau_min)"
)

CONVENTIONS = (
    "thiserror errors; no unwrap() in lib; rayon parallel (deterministic); "
    "f32 pipeline; feature='full-io' gates NetCDF; "
    "COG(≤50 bands) / GeoZarr-V3(>50 bands); "
    "band naming: rhot_443 / rhos_443 (nm suffix)"
)

TOLERANCES = (
    "ρt exact(1e-6) | ρs<0.002 | AOT<0.001 | Pearson R>0.999 | lat/lon exact(1e-8)"
)

COMMANDS = (
    "cargo test --features full-io | "
    "pytest tests/regression/ -v | "
    "cargo bench --bench performance"
)

# Dash-free versions for /compact output (Kiro rejects - and -- in the prompt)
PHYSICS_COMPACT = (
    "rhos=(rhot/tt_gas minus romix)/(dutott+astot*(rhot/tt_gas minus romix)) | "
    "DSF=intercept(200px darkest) | wave_range=(400,900)nm | tgas_min=0.85 | "
    "AOT=min RMSD model selection | NaN when ds<romix(tau_min)"
)

CONVENTIONS_COMPACT = (
    "thiserror errors; no unwrap() in lib; rayon parallel (deterministic); "
    "f32 pipeline; feature=full_io gates NetCDF; "
    "COG(50 bands or fewer) / GeoZarr V3(more than 50 bands); "
    "band naming: rhot_443 / rhos_443 (nm suffix)"
)

TOLERANCES_COMPACT = (
    "rho_t exact(1e-6) | rho_s<0.002 | AOT<0.001 | Pearson R>0.999 | lat/lon exact(1e-8)"
)

COMMANDS_COMPACT = (
    "cargo test (features=full_io) | "
    "pytest tests/regression/ | "
    "cargo bench"
)


def build_done_block(state: dict) -> str:
    loaders_done = [s for s, v in state["loaders"].items() if v]
    ac_done = [m for m, v in state["ac_modules"].items() if v]
    writers_done = [w for w, v in state["writers"].items() if v]
    extras = []
    if state["ancillary_wired"]:
        extras.append("ancillary-wired")
    if state["roi_wired"]:
        extras.append("roi-limit")
    if state["ancillary_interp_done"]:
        extras.append("ancillary-interp")
    if state["dem_pressure"]:
        extras.append("dem-pressure")
    if state["settings_file"]:
        extras.append("settings-file")
    if state["l2w_exists"]:
        extras.append("l2w-chl_oc3")
    if state["perf_ci"]:
        extras.append("perf-ci")
    return (
        f"LOADERS: {','.join(loaders_done) or 'none'} | "
        f"AC: {','.join(ac_done) or 'none'} | "
        f"WRITERS: {','.join(writers_done) or 'none'} | "
        f"EXTRAS: {','.join(extras) or 'none'} | "
        f"TESTS: {state['rust_tests']}Rust+{state['python_tests']}Py"
    )


def build_pending_block(state: dict) -> str:
    pending = [k for k in PRIORITY if not is_done(k, state)]
    lines = []
    for i, key in enumerate(pending[:6], 1):
        t = NEXT_TASKS[key]
        lines.append(f"NEXT[{i}] {key}: {t['short']}")
    return " | ".join(lines)


def build_compact_prompt(state: dict) -> str:
    """Build a /compact slash-command prompt preserving critical roadmap state.

    Paste the output directly into a Kiro chat session when the context
    window is getting full.  All dashes are removed from the preservation
    block because Kiro's /compact command rejects - and -- in the prompt.
    """
    loaders_done = [s for s, v in state["loaders"].items() if v]
    ac_done      = [m for m, v in state["ac_modules"].items() if v]
    writers_done = [w for w, v in state["writers"].items() if v]

    extras = []
    if state["ancillary_wired"]: extras.append("ancillary_wired")
    if state["roi_wired"]:       extras.append("roi_limit")
    if state["ancillary_interp_done"]: extras.append("ancillary_interp")
    if state["dem_pressure"]:    extras.append("dem_pressure")
    if state["settings_file"]:   extras.append("settings_file")
    if state["l2w_exists"]:      extras.append("l2w_chl_oc3")
    if state["perf_ci"]:         extras.append("perf_ci")

    # Replace any residual hyphens in sensor/module names with underscores
    def nodash(s: str) -> str:
        return s.replace("-", "_")

    done_line = (
        f"loaders({','.join(nodash(s) for s in loaders_done) or 'none'}) | "
        f"ac({','.join(nodash(m) for m in ac_done) or 'none'}) | "
        f"writers({','.join(nodash(w) for w in writers_done) or 'none'}) | "
        f"{','.join(extras) or 'no_extras'} | "
        f"{state['rust_tests']}Rust+{state['python_tests']}Py tests"
    )

    pending = [k for k in PRIORITY if not is_done(k, state)]
    next_lines = []
    for i, key in enumerate(pending[:5], 1):
        t = NEXT_TASKS[key]
        next_lines.append(f"NEXT[{i}] {nodash(key)}: {t['short']}.")

    lines = [
        "/compact Keep:",
        "ACOLITE-RS Rust port of Python ACOLITE atmospheric correction.",
        f"DONE: {done_line}",
        f"PHYSICS: {PHYSICS_COMPACT}",
        f"CONVENTIONS: {CONVENTIONS_COMPACT}",
        f"TOLERANCES: {TOLERANCES_COMPACT}",
    ] + next_lines + [
        f"VERIFY: cargo test (features=full_io) && pytest tests/regression/",
    ]
    return "\n".join(lines)


def build_context(state: dict, next_key: str | None = None) -> str:
    """Build the full compact context string for Kiro."""
    done_block    = build_done_block(state)
    pending_block = build_pending_block(state)

    header = (
        "=== ACOLITE-RS RESUME CONTEXT ===\n"
        "Rust port of Python ACOLITE atmospheric correction toolkit.\n"
        "Repo: acolite/ (Python ref) | src/ (Rust port) | tests/regression/ (pytest)\n"
        "Pipeline: loader → ac(calibration→gas→rayleigh→DSF) → writer\n"
    )

    body = (
        f"DONE: {done_block}\n"
        f"PHYSICS: {PHYSICS}\n"
        f"CONVENTIONS: {CONVENTIONS}\n"
        f"TOLERANCES: {TOLERANCES}\n"
        f"COMMANDS: {COMMANDS}\n"
        f"PENDING: {pending_block}\n"
    )

    if next_key:
        if next_key not in NEXT_TASKS:
            print(f"Unknown --next key '{next_key}'. Valid: {', '.join(NEXT_TASKS)}",
                  file=sys.stderr)
            sys.exit(1)
        task = NEXT_TASKS[next_key]
        task_block = (
            f"\n=== TASK ===\n"
            f"{task['detail']}\n"
            f"\n"
            f"After completing the code changes:\n"
            f"1. Run: cargo test --features full-io\n"
            f"2. Run: pytest tests/regression/ -v\n"
            f"3. Fix any failures before finishing.\n"
        )
    else:
        task_block = ""

    return header + body + task_block


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a compact Kiro resume context for the ACOLITE Rust port"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--next", metavar="TASK",
        choices=list(NEXT_TASKS.keys()),
        help=f"Include full detail for this next task ({', '.join(NEXT_TASKS.keys())})",
    )
    group.add_argument(
        "--auto", action="store_true",
        help="Automatically select the highest-priority incomplete task",
    )
    group.add_argument(
        "--compact", action="store_true",
        help="Emit a /compact slash-command prompt to paste into Kiro chat",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List all tasks with their current completion status and exit",
    )
    args = parser.parse_args()

    state = detect_state()

    if args.list:
        print("ACOLITE-RS task status (auto-detected from repo):\n")
        for key in PRIORITY:
            done = is_done(key, state)
            mark = "✅" if done else "⬜"
            print(f"  {mark} {key:15s} {NEXT_TASKS[key]['short']}")
        print()
        done_count = sum(1 for k in PRIORITY if is_done(k, state))
        print(f"  {done_count}/{len(PRIORITY)} phase-D/E/F tasks complete")
        print(f"  Tests: {state['rust_tests']} Rust, {state['python_tests']} Python")
        return

    if args.compact:
        print(build_compact_prompt(state))
        return

    next_key = args.next
    if args.auto:
        next_key = next((k for k in PRIORITY if not is_done(k, state)), None)
        if next_key is None:
            print("All tracked tasks complete.", file=sys.stderr)
            sys.exit(0)

    print(build_context(state, next_key))


if __name__ == "__main__":
    main()
