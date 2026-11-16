#!/usr/bin/env bash
# Paper 1 experiment sequence.  Run from code/embeddings/.
# Everything downloads on first use; models go to hf_cache, CIFAR to
# torch_cache, results to runs/.
set -euo pipefail

CONFIG="${1:-config.yaml}"
PY="${PYTHON:-python}"

mkdir -p runs hf_cache torch_cache

# 1. CLIP image embeddings on CIFAR-10 (image tower):
#    intrinsic dimension, curvature, RD curve, gain, outage.
$PY run_all.py --config "$CONFIG" --model clip_vit_b32 --dataset cifar10 \
    --stages extract,geometry,rd,gain,outage

# 2. CLIP text embeddings on the CIFAR-10 class prompts (text tower).
$PY run_all.py --config "$CONFIG" --model clip_vit_b32 --dataset cifar10 \
    --tag clip_vit_b32__cifar10_text --stages extract,geometry,rd

# 3. A second image encoder on the same data, for the mismatch experiment.
$PY run_all.py --config "$CONFIG" --model dinov2_s --dataset cifar10 \
    --stages extract,geometry,rd,gain

# 4. Mismatch between CLIP and DINOv2, matched by Procrustes.
$PY run_all.py --config "$CONFIG" --model clip_vit_b32 --dataset cifar10 \
    --stages mismatch

# 5. BERT sentence embeddings on AG News (pure text manifold).
$PY run_all.py --config "$CONFIG" --model minilm --dataset agnews \
    --stages extract,geometry,rd,gain,outage,task

# 6. Larger CLIP backbone, as a robustness check on D and curvature.
$PY run_all.py --config "$CONFIG" --model clip_vit_l14 --dataset cifar100 \
    --stages extract,geometry,rd

echo "done. reports under runs/*/report.json"
