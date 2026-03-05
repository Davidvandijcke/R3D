#!/bin/bash
#
# 02_mc_stata_parallel.sh — Launch parallel Stata MC workers
#
# Usage: bash 02_mc_stata_parallel.sh [NCORES]
#   NCORES: Number of parallel workers (default: 15)
#
# Each worker handles a chunk of simulations (1-100).
# Output files are uniquely named per sim, so no write conflicts.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Load Stata module
module load stata-mp/18 2>/dev/null || true

# Find Stata executable (try stata-mp first, fall back to stata-se, then stata)
if command -v stata-mp &>/dev/null; then
    STATA="stata-mp"
elif [ -x "/sw/pkgs/arc/stata-mp/18/stata-mp" ]; then
    STATA="/sw/pkgs/arc/stata-mp/18/stata-mp"
elif command -v stata-se &>/dev/null; then
    STATA="stata-se"
elif [ -x "/sw/pkgs/arc/stata-mp/18/stata-se" ]; then
    STATA="/sw/pkgs/arc/stata-mp/18/stata-se"
elif command -v stata &>/dev/null; then
    STATA="stata"
elif [ -x "/sw/pkgs/arc/stata-mp/18/stata" ]; then
    STATA="/sw/pkgs/arc/stata-mp/18/stata"
else
    echo "ERROR: Cannot find Stata executable"
    exit 1
fi

echo "Using Stata: $STATA"

SIMS=100
NCORES=${1:-15}
CHUNK=$((SIMS / NCORES))

echo "Launching $NCORES parallel workers for $SIMS simulations (chunk size ~$CHUNK)"
echo "Working directory: $SCRIPT_DIR"

# Create output directory
mkdir -p output/results_Stata

PIDS=()
for i in $(seq 0 $((NCORES - 1))); do
    START=$((i * CHUNK + 1))
    if [ $i -eq $((NCORES - 1)) ]; then
        END=$SIMS
    else
        END=$(((i + 1) * CHUNK))
    fi

    echo "  Worker $((i+1)): sims $START-$END"

    # Run Stata in batch mode with arguments
    # The -b flag runs in batch mode, arguments are passed after the .do file
    $STATA -b do 02_mc_stata_worker.do $START $END &
    PIDS+=($!)

    # Rename the default log file to avoid collisions
    # Stata -b creates 02_mc_stata_worker.log by default
    # But each worker also creates its own log inside output/results_Stata/
    # Small delay to avoid race condition on log file creation
    sleep 0.5
done

echo ""
echo "All $NCORES workers launched. PIDs: ${PIDS[*]}"
echo "Waiting for all workers to complete..."

# Wait for all workers and track failures
FAILED=0
for i in "${!PIDS[@]}"; do
    if ! wait "${PIDS[$i]}"; then
        echo "  WARNING: Worker $((i+1)) (PID ${PIDS[$i]}) failed"
        FAILED=$((FAILED + 1))
    fi
done

if [ $FAILED -gt 0 ]; then
    echo ""
    echo "WARNING: $FAILED worker(s) failed. Check logs in output/results_Stata/"
else
    echo ""
    echo "All $NCORES workers completed successfully!"
fi

# Count output files
N_RES=$(ls -1 output/results_Stata/mc_res_*.csv 2>/dev/null | wc -l)
N_PVAL=$(ls -1 output/results_Stata/mc_pval_*.csv 2>/dev/null | wc -l)
echo "Result files: $N_RES tau/cb CSVs, $N_PVAL p-value CSVs"
