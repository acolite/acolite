# ACOLITE Regression Test Harness
#
# Compares Python ACOLITE output against Rust acolite-rs output for PACE OCI.
# Ensures the Rust port produces numerically consistent results.
#
# Usage:
#   pytest tests/regression/ -v
#   pytest tests/regression/ -v -k pace
#   pytest tests/regression/ -v --runslow   # include tests requiring real data
