#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-${SCRIPT_DIR}/config.yaml}"
PY="${PYTHON:-python}"
"$PY" "$SCRIPT_DIR/run_all.py" --config "$CONFIG" --model clip_vit_b32 --dataset cifar10 --stages extract,geometry,rd,gain,outage
"$PY" "$SCRIPT_DIR/run_all.py" --config "$CONFIG" --model dinov2_s --dataset cifar10 --stages extract,geometry,rd,gain
"$PY" "$SCRIPT_DIR/run_all.py" --config "$CONFIG" --model clip_vit_b32 --dataset cifar10 --stages mismatch
"$PY" "$SCRIPT_DIR/run_all.py" --config "$CONFIG" --model minilm --dataset agnews --stages extract,geometry,rd,gain,outage,task
