#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data/dragon_stand}"
PREPARED_DIR="${PREPARED_DIR:-$SCRIPT_DIR/prepared}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
CONFIG="${CONFIG:-$SCRIPT_DIR/dragon_stand.ini}"
GMMREG_EXE="${GMMREG_EXE:-$REPO_ROOT/C++/build/gmmreg_demo}"

# ── Step 1: prepare point clouds and GT poses ─────────────────────────────────
echo "=== Step 1: prepare_data ==="
python3 "$SCRIPT_DIR/prepare_data.py" \
    --data_dir   "$DATA_DIR" \
    --output_dir "$PREPARED_DIR"

# ── Step 2: pairwise rigid registration and evaluation ────────────────────────
echo ""
echo "=== Step 2: run_experiment ==="
python3 "$SCRIPT_DIR/run_experiment.py" \
    --prepared_dir "$PREPARED_DIR" \
    --output_dir   "$RESULTS_DIR" \
    --config_tmpl  "$CONFIG" \
    --gmmreg_exe   "$GMMREG_EXE"
