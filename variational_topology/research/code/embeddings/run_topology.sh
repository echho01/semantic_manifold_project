#!/usr/bin/env bash
# Paper 2 experiment sequence.  Run from code/embeddings/ AFTER run_all.sh
# (which caches the embeddings).  CPU only.
set -euo pipefail
CONFIG="${1:-config.yaml}"
PY="${PYTHON:-python}"

# 1. CIFAR-10 images: CLIP vs DINOv2 (same items) -- floor, dictionary curve,
#    reparametrization on real data, cliffs with the 10 classes as concepts,
#    Reeb overhead of a probe score.
$PY run_topology.py --config "$CONFIG" --dataset cifar10 --modality image \
    --models clip_vit_b32,dinov2_s

# 2. AG News: a single text model -- reparametrization, cliffs (4 classes), Reeb.
$PY run_topology.py --config "$CONFIG" --dataset agnews --models minilm \
    --stages reparam,cliff,reeb

echo "done. reports under runs/topology__*/report.json"
