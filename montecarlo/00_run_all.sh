#!/bin/bash
#
# 00_run_all.sh — Master script for R3D Monte Carlo comparison
#
# Usage: bash 00_run_all.sh [phase]
#   phase 0: Re-run equivalence tests
#   phase 1: R Monte Carlo
#   phase 2: Stata Monte Carlo
#   phase 3: Compare results
#   (no arg): Run all phases

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

STATA="/Applications/Stata/StataSE.app/Contents/MacOS/StataSE"
STATA_R3D_DIR="$SCRIPT_DIR/../stata_r3d"

PHASE="${1:-all}"

# ============================================================================
# Phase 0: Equivalence tests
# ============================================================================
if [[ "$PHASE" == "0" || "$PHASE" == "all" ]]; then
    echo "=== Phase 0: Re-run equivalence tests ==="
    cd "$STATA_R3D_DIR/tests"
    echo "  Regenerating reference data..."
    Rscript generate_reference_data.R
    echo "  Running Stata equivalence tests..."
    "$STATA" -b do test_equivalence.do
    if grep -q "FAILED" test_equivalence.log 2>/dev/null; then
        echo "  WARNING: Some equivalence tests FAILED. Check test_equivalence.log"
    else
        echo "  All equivalence tests passed."
    fi
    cd "$SCRIPT_DIR"
    echo ""
fi

# ============================================================================
# Phase 1: R Monte Carlo
# ============================================================================
if [[ "$PHASE" == "1" || "$PHASE" == "all" ]]; then
    echo "=== Phase 1: R Monte Carlo ==="
    Rscript 01_mc_R.R
    echo ""
fi

# ============================================================================
# Phase 2: Stata Monte Carlo
# ============================================================================
if [[ "$PHASE" == "2" || "$PHASE" == "all" ]]; then
    echo "=== Phase 2: Stata Monte Carlo ==="
    if [[ ! -x "$STATA" ]]; then
        echo "  Stata not found at $STATA. Skipping Phase 2."
        echo "  Run manually: $STATA -b do 02_mc_stata.do"
    else
        "$STATA" -b do 02_mc_stata.do
        echo "  Stata MC complete. Check 02_mc_stata.log for details."
    fi
    echo ""
fi

# ============================================================================
# Phase 3: Compare
# ============================================================================
if [[ "$PHASE" == "3" || "$PHASE" == "all" ]]; then
    echo "=== Phase 3: Comparison ==="
    Rscript 03_mc_compare.R
    echo ""
fi

echo "=== Done ==="
