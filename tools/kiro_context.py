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
    python tools/kiro_context.py --next sentinel3
    python tools/kiro_context.py --next roi
    python tools/kiro_context.py --next ancillary
    python tools/kiro_context.py --next l2r-output
    python tools/kiro_context.py --next l2w
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
        "sentinel3": _file_exists("src/loader/sentinel3.rs"),  # stub sensor def only
        "emit":      _file_exists("src/loader/emit.rs"),
        "prisma":    _file_exists("src/loader/prisma.rs"),
    }

    # sentinel3 loader check: sensor def exists but loader is a placeholder
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
    }

    ancillary_wired = ac_modules["ancillary"] and bool(
        _grep_count(r"ancillary", "src/pipeline.rs")
    )

    writers = {
        "cog":      _file_exists("src/writer/cog.rs"),
        "geozarr":  _file_exists("src/writer/geozarr.rs"),
        "netcdf_l2r": _file_exists("src/writer/netcdf_l2r.rs"),
    }

    rust_tests   = _grep_count(r"#\[test\]", "src/**/*.rs", "tests/*.rs")
    python_tests = _grep_count(r"def test_", "tests/regression/**/*.py")

    # Pipeline feature flags from pipeline.rs
    roi_wired = bool(_grep_count(r"\blimit\b", "src/pipeline.rs"))
    l2w_exists = _file_exists("src/ac/l2w.rs")
    perf_ci    = _file_exists(".github/workflows/perf.yml")

    return {
        "loaders":          loaders,
        "ac_modules":       ac_modules,
        "ancillary_wired":  ancillary_wired,
        "writers":          writers,
        "roi_wired":        roi_wired,
        "l2w_exists":       l2w_exists,
        "perf_ci":          perf_ci,
        "rust_tests":       rust_tests,
        "python_tests":     python_tests,
    }


# ---------------------------------------------------------------------------
# Next-task definitions — ordered by roadmap priority
# ---------------------------------------------------------------------------

NEXT_TASKS = {
    "sentinel3": {
        "short": "Port Sentinel-3 OLCI NetCDF loader",
        "detail": (
            "Port acolite/sentinel3/l1_convert.py → src/loader/sentinel3.rs. "
            "Sensor def exists at src/sensors/sentinel3.rs (has band names + wavelengths, stub metadata). "
            "Follow src/loader/pace.rs for NetCDF reading (feature='full-io'). "
            "Parse xfdumanifest.xml for scene datetime/geometry; "
            "read Oa01_radiance.nc … Oa21_radiance.nc (21 bands); "
            "apply tie-point geometry from geo_coordinates.nc. "
            "Produce BandData arrays; hook into pipeline.rs sensor dispatch. "
            "Output: per-band COG (21 bands ≤ 50). "
            "Add tests/sentinel3_e2e.rs with a synthetic NetCDF fixture. "
            "Add tests/regression/test_sentinel3_regression.py (Tier 1 synthetic tests)."
        ),
    },
    "roi": {
        "short": "Implement --limit ROI subsetting for all loaders",
        "detail": (
            "Add limit=[S,W,N,E] parameter to all three loaders and to src/main.rs CLI. "
            "Landsat/S2 (GeoTIFF): use GDAL windowed reads — compute pixel window from "
            "bounding box and pass as GdalReadWindow to geotiff.rs read calls. "
            "PACE (NetCDF): slice latitude/longitude arrays to find row/col index range, "
            "then read only that subarray from each band variable. "
            "Add --limit S,W,N,E CLI arg parsed as four f64 values. "
            "Regression: compare Rust ROI-limited output vs Python ACOLITE with same limit."
        ),
    },
    "ancillary": {
        "short": "Wire ancillary data (ozone/pressure/wind) into pipeline",
        "detail": (
            "src/ac/ancillary.rs (180 lines) exists but is never called from src/pipeline.rs. "
            "Integrate: at AC start, call ancillary::fetch(datetime, lat, lon, credentials) "
            "using the scene centre coordinates from loader Metadata. "
            "Pass returned Ancillary{uoz, uwv, pressure, wind} into dsf::run() replacing "
            "hardcoded defaults (uoz=0.3, uwv=1.5, pressure=1013.25). "
            "Add --no-ancillary CLI flag to skip download and use defaults. "
            "Test: assert that a non-default ozone value changes rhos by a measurable amount."
        ),
    },
    "l2r-output": {
        "short": "Add NetCDF L2R output writer matching Python ACOLITE format",
        "detail": (
            "Create src/writer/netcdf_l2r.rs behind #[cfg(feature='full-io')]. "
            "Schema: global attrs (sensor, datetime, acolite_version, limit), "
            "dimensions (y, x), variables (lat f64, lon f64, rhos_NNN f32 per band, "
            "rhot_NNN f32 per band, wavelength f32). "
            "Match Python ACOLITE's variable naming exactly: rhos_443, rhot_443 etc. "
            "Add write_netcdf_l2r() to src/writer/mod.rs dispatch. "
            "Round-trip test: Python netCDF4 can open and read the Rust output."
        ),
    },
    "l2w": {
        "short": "Port first water product: chlorophyll-a OC3",
        "detail": (
            "Port OC3 chlorophyll-a from acolite/ac/parameters/chlorophyll.py → src/ac/l2w.rs. "
            "Function signature: fn chl_oc3(rhos: &HashMap<u32, Array2<f32>>) -> Array2<f32>. "
            "Band ratio: R = log10(max(rhos_443, rhos_490) / rhos_560). "
            "OC3 polynomial: chl = 10^(a0 + a1*R + a2*R^2 + a3*R^3 + a4*R^4) with "
            "NASA OC3 coefficients [0.2515, -2.3798, 1.5317, -0.1428, -0.9822]. "
            "Add --l2w-parameters chl_oc3 CLI flag. "
            "Write output as raster band in COG/GeoZarr with units='mg m-3'. "
            "Regression tests: known rhos values → expected Chl-a within 1e-4."
        ),
    },
    "perf-ci": {
        "short": "Add GitHub Actions performance regression gate",
        "detail": (
            "Create .github/workflows/perf.yml: trigger on push to main and PR. "
            "Steps: checkout → cargo bench --bench performance -- --output-format bencher | tee bench.txt. "
            "Store baselines in benches/performance_baseline.json: "
            "{\"landsat_mpx_s\": 940000, \"s2_mpx_s\": 576000, \"pace_mpx_s\": 18700000}. "
            "Parse bench.txt, compare to baseline, fail if any metric drops >10%. "
            "On main-branch merge: update baseline JSON and commit via bot. "
            "Also add memory check: /usr/bin/time -v cargo run … 2>&1 | grep 'Maximum resident'."
        ),
    },
}

# Priority order for --auto
PRIORITY = ["sentinel3", "roi", "ancillary", "l2r-output", "l2w", "perf-ci"]


def is_done(key: str, state: dict) -> bool:
    if key == "sentinel3":
        return state["loaders"].get("sentinel3", False)
    if key == "roi":
        return state["roi_wired"]
    if key == "ancillary":
        return state["ancillary_wired"]
    if key == "l2r-output":
        return state["writers"].get("netcdf_l2r", False)
    if key == "l2w":
        return state["l2w_exists"]
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


def build_done_block(state: dict) -> str:
    loaders_done = [s for s, v in state["loaders"].items() if v]
    ac_done = [m for m, v in state["ac_modules"].items() if v]
    writers_done = [w for w, v in state["writers"].items() if v]
    extras = []
    if state["ancillary_wired"]:
        extras.append("ancillary-wired")
    if state["roi_wired"]:
        extras.append("roi-limit")
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
    window is getting full.  The preservation instruction is ~400 tokens
    so it survives any compaction budget.
    """
    loaders_done = [s for s, v in state["loaders"].items() if v]
    ac_done      = [m for m, v in state["ac_modules"].items() if v]
    writers_done = [w for w, v in state["writers"].items() if v]

    extras = []
    if state["ancillary_wired"]: extras.append("ancillary-wired")
    if state["roi_wired"]:       extras.append("roi-limit")
    if state["l2w_exists"]:      extras.append("l2w-chl_oc3")
    if state["perf_ci"]:         extras.append("perf-ci")

    done_line = (
        f"loaders({','.join(loaders_done) or 'none'}) | "
        f"ac({','.join(ac_done) or 'none'}) | "
        f"writers({','.join(writers_done) or 'none'}) | "
        f"{','.join(extras) or 'no-extras'} | "
        f"{state['rust_tests']}Rust+{state['python_tests']}Py tests"
    )

    pending = [k for k in PRIORITY if not is_done(k, state)]
    next_lines = []
    for i, key in enumerate(pending[:5], 1):
        t = NEXT_TASKS[key]
        next_lines.append(f"NEXT[{i}] {key}: {t['short']}.")

    lines = [
        "/compact Keep:",
        f"ACOLITE-RS Rust port of Python ACOLITE atmospheric correction.",
        f"DONE: {done_line}",
        f"PHYSICS: {PHYSICS}",
        f"CONVENTIONS: {CONVENTIONS}",
        f"TOLERANCES: {TOLERANCES}",
    ] + next_lines + [
        f"VERIFY: cargo test --features full-io && pytest tests/regression/ -v",
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
