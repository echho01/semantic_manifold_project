# Embedding extraction and empirical diagnostics

This optional pipeline processes CLIP image/text features, MiniLM/BERT text
features and DINOv2 image features. The published synthetic calculations are
reproduced separately by `../synthetic/run_paper1.py` and `run_paper2.py`.

## Installation and execution

```bash
python -m pip install -r code/embeddings/requirements.txt
python code/embeddings/run_all.py --config code/embeddings/config.yaml \
  --model clip_vit_b32 --dataset cifar10 --stages extract,geometry,rd,gain,outage
```

Run another encoder on the identical dataset and then compare their paired
features:

```bash
python code/embeddings/run_all.py --config code/embeddings/config.yaml \
  --model dinov2_s --dataset cifar10 --stages extract,geometry,rd
python code/embeddings/run_all.py --config code/embeddings/config.yaml \
  --model clip_vit_b32 --dataset cifar10 --stages mismatch
```

Paths in the configuration are relative to the **configuration file**, so these
commands work from any current directory. The driver writes
`runs/<model>__<dataset>/arrays.npz` and `report.json`. Model extraction can use
CUDA, MPS or CPU. Subsequent stages use NumPy/scikit-learn and existing arrays.
Single-stage reruns preserve the other completed report sections.

## Measurement definitions

| Stage | What is measured |
|---|---|
| `extract` | Float features and stable original-item IDs; CLIP features are L2-normalized |
| `geometry` | TwoNN dimension and local Monge-chart estimates of `Scal`, `||II||²`, `||tr II||²` |
| `rd` | Held-out ambient MSE of fixed-rate k-means codebooks; rate is `log2(N)` bits |
| `gain` | Held-out noisy embeddings decoded with local affine planes or nearest reference samples |
| `outage` | Held-out nearest-centroid classification error, with clean error reported |
| `task` | Held-out quantization of a rank-specified supervised linear probe |
| `mismatch` | A-trained encoding indices decoded by B-space cluster centroids, compared with B-trained quantization |

The same deterministic 80/20 fit/evaluation split is used within each dataset.
Codebook sizes exceeding the training count are skipped and recorded without
changing their reported rate. Probes and centroids use training labels only.
Cross-encoder experiments require matching item IDs and ordering.

Fixed-rate k-means residuals are labelled as empirical diagnostics. A fitted
intercept cannot convert their second-order term into the Shannon RD
coefficient. The geometric Shannon coefficients are reported separately in
nats and refer to manifold-constrained reconstruction. Off-manifold k-means
centroids are evaluated using their actual ambient squared error.

The gain decoder receives the noisy observation and independent reference
samples. It projects onto a **centred** local plane with the estimated number of
tangent directions. Its zero-noise interpolation error is included. The ideal
`10 log10(D_ambient/d)` tangent reference is labelled separately.

Local curvature radii describe the fitted chart and are not treated as global
reach estimates. TwoNN removes exact duplicate rows and records their number.
CIFAR class prompts are saved once as `text_prototypes`; they are not repeated
into a pseudo-manifold of thousands of identical text vectors. STS sentence
pair scores are stored as `pair_scores` with `pair_indices`, separate from
sample-level classification labels.

## CPU verification

```bash
python code/embeddings/selftest.py
python code/embeddings/selftest_rd.py
```

These checks exercise high-ambient-dimensional geometry, translation-invariant
projection, probe rank, paired-edge metrics, quantizer transfer, continuous RD
coefficients, JSON serialization, staged resumption and standalone outage runs.
No pretrained weights are downloaded by the tests.
