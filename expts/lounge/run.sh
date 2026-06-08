#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data/lounge}"
PREPARED_DIR="${PREPARED_DIR:-$SCRIPT_DIR/prepared}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
CONFIG="${CONFIG:-$REPO_ROOT/expts/lounge.ini}"
SEQ_RIGID3D="${SEQ_RIGID3D:-$REPO_ROOT/C++/build/seq_rigid3d}"

# ── Step 1: prepare point clouds and GT poses ─────────────────────────────────
echo "=== Step 1: prepare_data ==="
python3 "$SCRIPT_DIR/prepare_data.py" \
    --data_dir   "$DATA_DIR" \
    --output_dir "$PREPARED_DIR"

# ── Step 2: sequential pairwise rigid registration ────────────────────────────
echo ""
echo "=== Step 2: seq_rigid3d ==="
"$SEQ_RIGID3D" \
    --pcloud_dir "$PREPARED_DIR/pclouds" \
    --output_dir "$RESULTS_DIR" \
    --config     "$CONFIG" \
    --gt_dir     "$PREPARED_DIR/gt_poses" \
    --step       5 \
    --max_pairs  -1

# ── Results ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Results ==="
cat "$RESULTS_DIR/results.json"
