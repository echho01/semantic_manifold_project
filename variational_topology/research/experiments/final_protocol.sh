#!/bin/bash
# Final paired evaluation: three arms fine-tuned from one converged MSE codec per seed.
cd /mnt/workspace/gaofan/projects/semantic-manifold-topology
export PYTHONPATH=$PWD/pydeps OMP_NUM_THREADS=4
mkdir -p runs/final
T="python3 experiments/topology_jscc.py"
COMMON="--device cuda --data data/omniglot --batch 128 --size 32 --rho 0.4 --groups 4 --workers 8 --hidden 64"
EVAL="--device cuda --data data/omniglot --split test --limit 2048 --seed 1115 --snr 0 10 20 --channels awgn rayleigh --repeats 1 --queries 64 --batch 64"
one() {  # seed gpu
  s=$1; export CUDA_VISIBLE_DEVICES=$2
  $T train $COMMON --seed $s --mode mse --epochs 100 --lr 3e-4 --checkpoint runs/final/base_$s.pt || return 1
  FT="$COMMON --seed $s --epochs 20 --lr 1e-4 --warmup 0 --init runs/final/base_$s.pt"
  ( $T train $FT --mode mse --checkpoint runs/final/mse_$s.pt && $T test $EVAL --checkpoint runs/final/mse_$s.pt --output runs/final/mse_$s.json ) &
  ( $T train $FT --mode query --lambda-query 0.01 --structure-images 32 --queries 128 --checkpoint runs/final/query_$s.pt && $T test $EVAL --checkpoint runs/final/query_$s.pt --output runs/final/query_$s.json ) &
  ( $T train $FT --mode proposed --lambda-query 0.01 --lambda-mst 0.02 --structure-images 32 --queries 128 --checkpoint runs/final/proposed_$s.pt && $T test $EVAL --checkpoint runs/final/proposed_$s.pt --output runs/final/proposed_$s.json ) &
  ( $T train $FT --mode topojscc --structure-images 32 --checkpoint runs/final/topojscc_$s.pt && $T test $EVAL --checkpoint runs/final/topojscc_$s.pt --output runs/final/topojscc_$s.json ) &
  wait
}
one 1115 0 > logs/final_1115.log 2>&1 &
one 1116 1 > logs/final_1116.log 2>&1 &
one 1117 0 > logs/final_1117.log 2>&1 &
wait
echo ALLDONE > logs/final.done
