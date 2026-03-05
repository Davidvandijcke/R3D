#!/bin/bash
#
# run_stata_mc_nohup.sh — Launch parallel Stata MC workers under nohup
#
# Usage: nohup bash run_stata_mc_nohup.sh &
# Or just: bash run_stata_mc_nohup.sh  (it wraps itself in nohup-style)
#
# Output goes to a timestamped directory under montecarlo/runs/

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Timestamp for this run
TS=$(date +%Y%m%d_%H%M%S)
RUNDIR="$SCRIPT_DIR/runs/${TS}"
mkdir -p "$RUNDIR"

# Master log
LOGFILE="$RUNDIR/master.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== Stata MC Parallel Run ==="
echo "Timestamp: $TS"
echo "Run directory: $RUNDIR"
echo "Started: $(date)"
echo ""

# Load Stata module
source /etc/profile.d/lmod.sh 2>/dev/null || true
module load stata-mp/18 2>/dev/null || true

# Find Stata executable
for candidate in stata-mp stata-se stata; do
    if command -v $candidate &>/dev/null; then
        STATA=$candidate
        break
    fi
done
if [ -z "$STATA" ]; then
    for candidate in /sw/pkgs/arc/stata-mp/18/stata-mp /sw/pkgs/arc/stata-mp/18/stata-se /sw/pkgs/arc/stata-mp/18/stata; do
        if [ -x "$candidate" ]; then
            STATA=$candidate
            break
        fi
    done
fi
if [ -z "$STATA" ]; then
    echo "ERROR: Cannot find Stata executable"
    exit 1
fi
echo "Using Stata: $STATA"

# Make sure output dirs exist (shared with R output)
mkdir -p output/results_Stata

# Symlink run dir for easy access to logs
ln -sfn "$RUNDIR" "$SCRIPT_DIR/runs/latest"

SIMS=100
NCORES=15
CHUNK=$((SIMS / NCORES))

echo "Launching $NCORES parallel workers for $SIMS simulations (chunk size ~$CHUNK)"
echo ""

PIDS=()
for i in $(seq 0 $((NCORES - 1))); do
    START=$((i * CHUNK + 1))
    if [ $i -eq $((NCORES - 1)) ]; then
        END=$SIMS
    else
        END=$(((i + 1) * CHUNK))
    fi

    WORKER_LOG="$RUNDIR/worker_${START}_${END}.log"
    echo "  Worker $((i+1)): sims $START-$END -> $WORKER_LOG"

    # Run Stata in batch mode; redirect its output to per-worker log
    $STATA -b do 02_mc_stata_worker.do $START $END > "$WORKER_LOG" 2>&1 &
    PIDS+=($!)
done

echo ""
echo "All $NCORES workers launched. PIDs: ${PIDS[*]}"
echo "PID file: $RUNDIR/pids.txt"
printf '%s\n' "${PIDS[@]}" > "$RUNDIR/pids.txt"
echo ""
echo "Waiting for all workers to complete..."
echo "(Safe to disconnect SSH - processes will continue)"
echo ""

# Wait for all workers and track failures
FAILED=0
for i in "${!PIDS[@]}"; do
    PID=${PIDS[$i]}
    if wait "$PID"; then
        echo "  Worker $((i+1)) (PID $PID) completed OK"
    else
        echo "  WARNING: Worker $((i+1)) (PID $PID) FAILED (exit code $?)"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Summary ==="
echo "Finished: $(date)"

if [ $FAILED -gt 0 ]; then
    echo "WARNING: $FAILED worker(s) failed. Check logs in $RUNDIR/"
else
    echo "All $NCORES workers completed successfully!"
fi

# Count output files
N_RES=$(ls -1 output/results_Stata/mc_res_*.csv 2>/dev/null | wc -l)
N_PVAL=$(ls -1 output/results_Stata/mc_pval_*.csv 2>/dev/null | wc -l)
echo "Result files: $N_RES tau/cb CSVs, $N_PVAL p-value CSVs"
echo ""
echo "Run directory: $RUNDIR"
echo "To check coverage, run: Rscript 03_mc_compare.R"
