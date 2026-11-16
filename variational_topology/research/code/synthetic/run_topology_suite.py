#!/usr/bin/env python
"""Reproduce Paper 2: measured topology, channel events, codes, and task entropy.

All radii use half the maximum edge length. Noise has TOTAL power sigma^2,
with per-coordinate variance sigma^2 / D_a. Every Monte Carlo replicate
uses one standard-noise cloud throughout its sigma sweep. Results contain
per-replicate rows, input identities, exact configurations, source hashes,
and dependency versions. The four vector figures are generated from these
saved measurements; no manuscript number is fitted by hand.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ncx2, norm
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path[:0] = [str(ROOT / 'code' / 'topology'), str(ROOT / 'code' / 'embeddings')]
import tda
import reeb
import sync

SEED = 20261020


def _json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Unsupported result value {type(value).__name__}")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=_json_value) + '\n')


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def wilson(success, n, z=1.959963984540054):
    p = success / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return float(center - half), float(center + half)


def distance(A, B):
    return max(tda.bottleneck_by_degree(A, B))


def _paired_rips(D):
    dg = tda.rips(D=D, maxdim=1)
    dg[0] = tda.h0_diagram(tda.h0_deaths(D=D))
    return dg


def dictionary(out):
    rows = []
    for kind, counts in [('two_circles', [4, 8, 12, 16, 24, 32, 36, 37, 38, 40, 48, 64]),
                         ('clifford_torus', [6, 10, 14, 18, 22, 24, 26, 27, 28, 30])]:
        reach = 0.5 if kind == 'two_circles' else 1.0
        truth = (2, 2) if kind == 'two_circles' else (1, 2, 1)
        eta = reach / 6
        b = (2 + np.sqrt(2)) * eta
        for m in counts:
            t = 2 * np.pi * np.arange(m) / m
            if kind == 'two_circles':
                X = np.concatenate([np.c_[s * 1.5 + np.cos(t), np.sin(t)] for s in (-1, 1)])
                eps = 2 * np.sin(np.pi / (2 * m))
            else:
                u, v = np.meshgrid(t, t)
                X = np.c_[np.cos(u.ravel()), np.sin(u.ravel()), np.cos(v.ravel()), np.sin(v.ravel())]
                eps = 2 * np.sqrt(2) * np.sin(np.pi / (2 * m))
            dg = tda.rips(X=X, maxdim=len(truth) - 1, thresh=1.01 * b)
            est = tuple(tda.rank_estimator(dg[k], eta, b) for k in range(len(truth)))
            rows.append(dict(family=kind, m=m, atoms=len(X), fixed_bits=math.ceil(math.log2(len(X))),
                cover_radius=float(eps), reach=reach, eta=eta, second_radius=float(b),
                beta0=est[0], beta1=est[1], beta2=est[2] if len(est) > 2 else '',
                correct=est == truth, certified=bool(eps < eta)))
    write_csv(out / 'dictionary.csv', rows)
    summary = {'rows': rows, 'certified_cases': sum(r['certified'] for r in rows),
        'certificate_violations': sum(r['certified'] and not r['correct'] for r in rows),
        'dictionary_class': 'regular constructions; counts are not optimal covering numbers'}
    assert summary['certificate_violations'] == 0
    return summary


def reparameterization(out, trials):
    n, harmonic, delta = 80, 3, 0.05
    t = 2 * np.pi * (np.arange(n) + 0.137) / n
    Z = np.concatenate([np.c_[np.cos(k * t), np.sin(k * t)] for k in range(1, harmonic + 1)], axis=1) / np.sqrt(harmonic)
    cleanD = squareform(pdist(Z))
    dX = _paired_rips(cleanD)
    references = {sigma: (tda.noise_reference(cleanD, sigma),
                  _paired_rips(tda.noise_reference(cleanD, sigma))) for sigma in [0.5, 1.0, 2.0]}
    rows, calibration = [], []
    for Da in [64, 256, 1024, 4096]:
        frame = tda.random_frame(Z.shape[1], Da, seed=SEED + Da)
        X = Z @ frame
        for trial in range(trials):
            seed = SEED + 100000 + Da * 100 + trial
            G = np.random.default_rng(seed).normal(size=X.shape) / np.sqrt(Da)
            for sigma in [0.5, 1.0, 2.0]:
                Y = X + sigma * G
                DY = np.sqrt(tda.sq_dists(Y))
                refD, dP = references[sigma]
                dY = _paired_rips(DY)
                dD = _paired_rips(tda.debias(Y, sigma))
                raw, deb = distance(dY, dP), distance(dD, dX)
                rows.append(dict(D_a=Da, sigma=sigma, trial=trial, seed=seed,
                    naive=distance(dY, dX), raw_error=raw, debiased_error=deb,
                    raw_bound=tda.zeta_raw(sigma, n, Da, delta),
                    debiased_bound=tda.zeta_debiased(sigma, n, Da, delta),
                    raw_edge_error=float(np.max(np.abs(DY - refD)) / 2),
                    max_noise=float(np.max(np.linalg.norm(sigma * G, axis=1)))))
                if Da == 1024 and sigma == 1:
                    a = np.sqrt(np.log(2 * n * n / delta) / Da)
                    alpha, beta = 4 * sigma * a, 4 * sigma * sigma * (a + a * a)
                    for relative_variance_error in [-0.2, -0.1, 0, 0.1, 0.2]:
                        estimated_variance = sigma * sigma * (1 + relative_variance_error)
                        dd = _paired_rips(tda.debias(Y, np.sqrt(estimated_variance)))
                        v = abs(estimated_variance - sigma * sigma)
                        bound = (alpha + np.sqrt(alpha * alpha + 4 * (beta + 2 * v))) / 4
                        calibration.append(dict(trial=trial, seed=seed,
                            relative_variance_error=relative_variance_error,
                            error=distance(dd, dX), bound=float(bound)))
    np.savez_compressed(out / 'reparam_clean.npz', intrinsic_cloud=Z)
    write_csv(out / 'reparameterization.csv', rows)
    write_csv(out / 'calibration.csv', calibration)
    summary = {'n': n, 'harmonics': harmonic, 'delta_per_cloud': delta, 'trials': trials,
        'rows': rows, 'calibration': calibration,
        'raw_violations': sum(r['raw_error'] > r['raw_bound'] + 2e-7 for r in rows),
        'debiased_violations': sum(r['debiased_error'] > r['debiased_bound'] + 2e-7 for r in rows),
        'calibration_violations': sum(r['error'] > r['bound'] + 2e-7 for r in calibration)}
    assert summary['raw_violations'] + summary['debiased_violations'] + summary['calibration_violations'] == 0
    return summary


def channel(out, trials):
    m, gap, theta = 64, 1.0, 0.25
    n = 2 * m
    rows, contact_rows, onsets = [], [], []
    dims = [64, 256, 1024, 4096]
    for Da in dims:
        scale = gap * (Da / np.log(m)) ** 0.25
        normalized_grid = np.unique(np.round(np.r_[np.linspace(0.15, 0.75, 17),
                                                    np.linspace(0.25, 0.45, 41)], 12))
        sigmas = np.unique(np.r_[0.0, normalized_grid * scale,
                                [0.25, np.sqrt(2) * theta, 0.5], np.linspace(0.30, 0.45, 16)])
        centers = np.r_[np.zeros(m), np.ones(m) * gap]
        center_diff = centers[:, None] - centers[None, :]
        baseline = center_diff ** 2
        replicate = []
        for trial in range(trials):
            seed = SEED + 200000 + Da * 1000 + trial
            G = np.random.default_rng(seed).normal(size=(n, Da)) / np.sqrt(Da)
            gramD2 = tda.sq_dists(G)
            cross_term = 2 * center_diff * (G[:, 0, None] - G[None, :, 0])
            first = {key: None for key in ['raw_spurious', 'zero_cross_contact', 'partition_failure', 'debiased_collapse']}
            for sigma in sigmas:
                D2 = np.maximum(baseline + sigma * cross_term + sigma * sigma * gramD2, 0)
                np.fill_diagonal(D2, 0)
                debD2 = np.maximum(D2 - 2 * sigma * sigma, 0)
                np.fill_diagonal(debD2, 0)
                raw_death = tda.h0_deaths(D=np.sqrt(D2))
                debD = np.sqrt(debD2)
                deb_death = tda.h0_deaths(D=debD)
                r_in = max(tda.h0_deaths(D=debD[:m, :m])[0], tda.h0_deaths(D=debD[m:, m:])[0])
                r_out = float(debD[:m, m:].min()) / 2
                raw_count = int(np.sum(raw_death > theta)) + 1
                deb_count = int(np.sum(deb_death > theta)) + 1
                true_partition = bool(r_in <= theta < r_out)
                zero_contact = bool(np.any(debD2[:m, m:] == 0))
                kb_predictions = centers + sigma * G[:, 0] > gap / 2
                accuracy = float(np.mean(kb_predictions == (centers > 0)))
                row = dict(D_a=Da, sigma=float(sigma), normalized_sigma=float(sigma / scale),
                    trial=trial, seed=seed, raw_count=raw_count, debiased_count=deb_count,
                    r_in=float(r_in), r_out=r_out, true_partition=true_partition,
                    zero_cross_contact=zero_contact, debiased_collapse=bool(deb_death[0] == 0),
                    raw_spurious=raw_count > 2, kb_accuracy=accuracy)
                rows.append(row); replicate.append(row)
                for name, event in [('raw_spurious', raw_count > 2), ('zero_cross_contact', zero_contact),
                                    ('partition_failure', not true_partition), ('debiased_collapse', deb_death[0] == 0)]:
                    if event and first[name] is None:
                        first[name] = float(sigma)
            onsets.append(dict(D_a=Da, trial=trial, seed=seed, **first))
        for sigma in sigmas:
            selected = [r for r in replicate if r['sigma'] == sigma]
            normalized = float(sigma / scale)
            if sigma == 0:
                pair_zero = pair_theta = 0.0
            else:
                noncentrality = Da * gap * gap / (2 * sigma * sigma)
                pair_zero = float(ncx2.cdf(Da, Da, noncentrality))
                pair_theta = float(ncx2.cdf(Da * (1 + 2 * theta * theta / (sigma * sigma)), Da, noncentrality))
            for event, getter, pair_p in [
                ('zero_cross_contact', lambda r: r['zero_cross_contact'], pair_zero),
                ('theta_cross_contact', lambda r: r['r_out'] <= theta, pair_theta),
                ('partition_failure', lambda r: not r['true_partition'], None),
                ('raw_spurious', lambda r: r['raw_spurious'], None),
                ('debiased_collapse', lambda r: r['debiased_collapse'], None)]:
                count = sum(getter(r) for r in selected)
                lo, hi = wilson(count, trials)
                lower = float(-np.expm1(m * np.log1p(-pair_p))) if pair_p is not None and pair_p < 1 else (1.0 if pair_p == 1 else None)
                upper = float(min(1, m * m * pair_p)) if pair_p is not None else None
                contact_rows.append(dict(D_a=Da, sigma=float(sigma), normalized_sigma=normalized,
                    event=event, successes=count, trials=trials, probability=count / trials,
                    ci_low=lo, ci_high=hi, pair_probability=pair_p, lower_bound=lower, upper_bound=upper))
    write_csv(out / 'channel_raw.csv', rows)
    write_csv(out / 'channel_probabilities.csv', contact_rows)
    write_csv(out / 'channel_onsets.csv', onsets)
    onset_summary = []
    for Da in dims:
        r = [x for x in onsets if x['D_a'] == Da]
        scale = gap * (Da / np.log(m)) ** 0.25
        summary = {'D_a': Da, 'noise_scale': float(scale),
                   'kb_per_point_one_percent_sigma': float(gap * np.sqrt(Da) / (2 * norm.isf(0.01)))}
        for name in ['raw_spurious', 'zero_cross_contact', 'partition_failure', 'debiased_collapse']:
            observed = [x[name] for x in r if x[name] is not None]
            values = observed + [np.inf] * (trials - len(observed))
            def quantile(q):
                value = np.quantile(values, q, method='inverted_cdf')
                return float(value) if np.isfinite(value) else None
            summary[name] = {'observed': len(observed), 'trials': trials,
                'right_censored': trials - len(observed), 'median': quantile(.5),
                'q025': quantile(.025), 'q975': quantile(.975)}
        onset_summary.append(summary)
    block_target = 0.10
    per_point_target = float(-np.expm1(np.log1p(-block_target) / n))
    block_comparison = []
    for Da in dims:
        probabilities = [x for x in contact_rows
            if x['D_a'] == Da and x['event'] == 'partition_failure']
        crossing_index = next((i for i, x in enumerate(probabilities)
            if x['probability'] >= block_target), None)
        crossing = probabilities[crossing_index] if crossing_index is not None else None
        below = probabilities[crossing_index - 1] if crossing_index is not None and crossing_index > 0 else None
        block_comparison.append({'D_a': Da, 'block_error_target': block_target,
            'monte_carlo_trials': trials,
            'graph_last_below_sigma': below['sigma'] if below else None,
            'graph_first_at_or_above_sigma': crossing['sigma'] if crossing else None,
            'graph_last_below_normalized_sigma': below['normalized_sigma'] if below else None,
            'graph_first_at_or_above_normalized_sigma': crossing['normalized_sigma'] if crossing else None,
            'graph_last_below_failure_probability': below['probability'] if below else None,
            'graph_first_at_or_above_failure_probability': crossing['probability'] if crossing else None,
            'graph_last_below_wilson95_low': below['ci_low'] if below else None,
            'graph_last_below_wilson95_high': below['ci_high'] if below else None,
            'graph_first_at_or_above_wilson95_low': crossing['ci_low'] if crossing else None,
            'graph_first_at_or_above_wilson95_high': crossing['ci_high'] if crossing else None,
            'known_center_block_sigma': float(gap * np.sqrt(Da) / (2 * norm.isf(per_point_target))),
            'known_center_per_point_error_at_threshold': per_point_target})
    write_csv(out / 'channel_block_reliability.csv', block_comparison)
    return {'points_per_concept': m, 'gap': gap, 'theta': theta, 'trials': trials,
        'block_reliability': block_comparison,
        'block_crossing_estimator': 'adjacent sampled sigmas bracketing the first marginal partition-failure frequency reaching 0.10; this is a grid bracket, not a confidence interval for the population threshold or a first-event onset quantile',
        'sigma_grid_definition': 'union of 17 normalized levels 0.15..0.75, 41 normalized levels 0.25..0.45 (step 0.005), 16 absolute levels 0.30..0.45, and 0,0.25,sqrt(2)/4,0.5',
        'onset_estimator': 'first sampled sigma per fixed noise path; empirical inverse-CDF quantiles include right-censored replicates; null quantile means beyond sweep',
        'onsets': onset_summary, 'probabilities': contact_rows,
        'event_definitions': {'zero_cross_contact': 'at least one debiased cross distance is zero',
            'debiased_collapse': 'the entire graph is connected at radius zero',
            'partition_failure': 'the theta graph does not recover both true labelled components',
            'raw_spurious': 'raw H0 threshold count exceeds two'}}


def mismatch_sync(out):
    n = 720
    t = 2 * np.pi * (np.arange(n) + 0.173) / n
    XT = np.c_[np.cos(t), np.sin(t)]
    s = np.where(t < np.pi, -1, 1)
    u = 2 * np.mod(t, np.pi)
    Xsplit = np.c_[s * 1.5 + np.cos(u), np.sin(u)]
    Xaccent = np.c_[np.cos(t), 0.9 * np.sin(t)]
    gT = tda.alpha2d(XT)
    rows = []
    order, _ = tda.maxmin(XT, 512, seed=SEED)
    for name, XR in [('split_components', Xsplit), ('ellipse_deformation', Xaccent)]:
        gR = tda.alpha2d(XR)
        delta = distance(gT, gR)
        for K in [8, 16, 32, 64, 128, 256, 512]:
            idx = order[:K]
            eps = float(cKDTree(XR[idx]).query(XR)[0].max())
            distortion = distance(gT, tda.alpha2d(XR[idx]))
            rows.append(dict(family=name, atoms=K, bits=math.ceil(math.log2(K)),
                reference_delta=delta, distortion=distortion, cover_radius=eps,
                deviation=abs(distortion - delta), within_band=abs(distortion - delta) <= eps + 1e-8))
    deformation = []
    for eps in [.4, .2, .1, .05, .02, .01]:
        XR = np.c_[np.cos(t), (1 - eps) * np.sin(t)]
        deformation.append({'epsilon': eps, 'diagram_distance': distance(gT, tda.alpha2d(XR))})
    interval_distance = distance(gT, tda.alpha2d(np.c_[np.cos(t), np.zeros(n)]))
    code_rows = []
    rng = np.random.default_rng(SEED + 300000)
    for K in [1, 2, 6, 12]:
        reference = np.c_[1 + np.arange(K), 2 + np.arange(K)]
        eps0 = .08
        for ratio in [2.0, 2.5, 4.5, 6.5, 9.5]:
            eps1 = eps0 / ratio
            for trial in range(25):
                offsets = rng.uniform(-eps0, eps0, size=(K, 2))
                if trial == 0:
                    offsets[:] = eps0
                if trial == 1:
                    offsets[:] = -eps0
                message, m = sync.encode(reference, reference + offsets, eps0, eps1)
                recovered = sync.decode(reference, message, eps0, m)
                error = tda.bottleneck(reference + offsets, recovered)
                expected = (m ** (2 * K) - 1).bit_length()
                code_rows.append(dict(K=K, ratio=ratio, m=m, trial=trial,
                    eps0=eps0, eps1=eps1, payload_bits=len(message), optimal_fixed_bits=expected,
                    coordinate_wise_bits=2 * K * math.ceil(math.log2(m)),
                    bottleneck_error=error, reconstruction_ok=error <= eps1 + 1e-12,
                    length_ok=len(message) == expected))
    write_csv(out / 'mismatch.csv', rows)
    write_csv(out / 'deformation.csv', deformation)
    write_csv(out / 'sync_codes.csv', code_rows)
    np.savez_compressed(out / 'mismatch_clouds.npz', transmitter=XT, split_receiver=Xsplit, accent_receiver=Xaccent, item_ids=np.arange(n))
    assert all(r['within_band'] for r in rows)
    assert all(r['reconstruction_ok'] and r['length_ok'] for r in code_rows)
    return {'n_reference_items': n, 'rows': rows, 'deformation': deformation,
        'true_interval_distance': interval_distance,
        'finite_below_limit_cases': [r for r in rows if r['distortion'] < r['reference_delta'] - 1e-8],
        'code_checks': len(code_rows), 'code_violations': 0, 'code_rows': code_rows,
        'sync_source': 'labelled separated local endpoint chart, births >= eps0; K/reference are shared'}


def torus_sample(n, seed):
    rng = np.random.default_rng(seed)
    u = rng.uniform(0, 2 * np.pi, n)
    v = []
    total = 0
    while total < n:
        candidate = rng.uniform(0, 2 * np.pi, n)
        candidate = candidate[rng.uniform(0, 3, n) < 2 + np.cos(candidate)]
        v.append(candidate); total += len(candidate)
    v = np.concatenate(v)[:n]
    radius = 2 + np.cos(v)
    return np.c_[radius * np.cos(u), radius * np.sin(u), np.sin(v)]


def entropy(labels):
    _, counts = np.unique(labels, return_counts=True, axis=0)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def task(out):
    # Integrate out the azimuth exactly before numerical quadrature.
    integrand = lambda v: (2 + np.cos(v)) / (4 * np.pi) * (2 / np.pi) * np.arcsin(1 / (2 + np.cos(v)))
    value, error = quad(integrand, 0, 2 * np.pi, points=[np.pi], epsabs=1e-12, epsrel=1e-12)
    X = torus_sample(500000, SEED + 400000)
    f = X[:, 0]
    middle = np.abs(f) < 1
    edge = np.where(middle, np.where(X[:, 1] > 0, 1, -1), np.where(f > 0, 2, -2))
    rows = []
    for step in [.4, .2, .1, .05, .02, .01]:
        q = np.floor(f / step).astype(np.int64)
        value_entropy = entropy(q)
        branch_entropy = entropy(np.c_[q, edge])
        rows.append(dict(step=step, measured_mse=float(np.mean((f - (q + .5) * step) ** 2)),
            score_entropy=value_entropy, branch_entropy=branch_entropy,
            excess=branch_entropy - value_entropy, exact_limit=float(value)))
    Xs = torus_sample(20000, SEED + 400001)
    diagnostics = []
    for link in [.25, .30, .35]:
        for n_bins in [20, 30, 40]:
            est = reeb.branch_entropy(Xs, Xs[:, 0], n_bins=n_bins, link=link, min_size=1)
            mp = reeb.mapper_graph(Xs, Xs[:, 0], n_int=12, overlap=.3, link=link, min_size=5)
            diagnostics.append(dict(link=link, n_bins=n_bins, binned_entropy=est,
                mapper_nodes=mp['nodes'], mapper_edges=mp['edges'], mapper_beta0=mp['beta0'], mapper_beta1=mp['beta1']))
    count = int(middle.sum())
    interval = wilson(count, len(X))
    write_csv(out / 'reeb_quantization.csv', rows)
    write_csv(out / 'reeb_diagnostics.csv', diagnostics)
    np.savez_compressed(out / 'reeb_cloud.npz', point_cloud=Xs)
    return {'source_measure': 'uniform surface area on torus of revolution, major radius 2 and minor radius 1',
        'n_entropy_samples': len(X), 'entropy_seed': SEED + 400000,
        'exact_branch_entropy': float(value), 'quadrature_error_estimate': float(error),
        'empirical_middle_mass': count / len(X), 'middle_mass_wilson95': interval,
        'quantization': rows, 'diagnostics': diagnostics}


def digits(out):
    from sklearn.datasets import load_digits
    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split
    from run_topology import run
    X, y = load_digits(return_X_y=True)
    X = X.astype(np.float64) / 16
    ids = np.asarray([f'digits:{i:04d}' for i in range(len(X))])
    reference, evaluation = train_test_split(np.arange(len(X)), test_size=.4, stratify=y, random_state=SEED)
    projection = PCA(n_components=16, svd_solver='full').fit(X[reference])
    Y = projection.transform(X)
    cache = out / 'digits_cache'
    for model, cloud in [('pixels64', X), ('pca16', Y)]:
        path = cache / f'{model}__digits'
        path.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path / 'arrays.npz', image=cloud.astype(np.float64), labels=y, item_ids=ids)
    cfg = {'seed': SEED, 'paths': {'root': str(cache)}, 'topology': {
        'n_items': 200, 'maxdim': 1, 'test_fraction': .4,
        'dictionary_sizes': [16, 32, 64, 128, 200], 'reparam_points': 80,
        'noise_levels': [.25, .5, 1.0], 'cliff_noise_levels': [0, .25, .5, 1.0, 2.0],
        'cliff_per_class': 20, 'reeb_bins': 12}}
    report = run(cfg, 'digits', ['pixels64', 'pca16'], ['floor', 'dictionary', 'reparam', 'cliff', 'reeb'])
    report['dataset'] = {'name': 'scikit-learn optical recognition of handwritten digits',
        'source': 'https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_digits.html',
        'n_items': len(X), 'reference_n': len(reference), 'evaluation_n': len(evaluation),
        'representation': '64 raw pixels and 16-dimensional PCA fitted only to reference IDs',
        'pretrained_models': False, 'pca_explained_variance_ratio': float(projection.explained_variance_ratio_.sum())}
    assert set(report['provenance']['reference_ids']) == set(ids[reference])
    assert set(report['provenance']['reference_ids']).isdisjoint(report['provenance']['evaluation_ids'])
    write_json(out / 'digits_report.json', report)
    return report


def figures(out, result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8,
        'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.top': False,
        'axes.spines.right': False, 'savefig.dpi': 180})
    palette = ['#2563a6', '#bd4c29', '#3b8054', '#855da6']

    def save(fig, name):
        fig.savefig(out / (name + '.pdf'), bbox_inches='tight')
        fig.savefig(out / (name + '.png'), bbox_inches='tight', dpi=180)
        plt.close(fig)

    if 'dictionary' in result and 'mismatch' in result:
        fig, axes = plt.subplots(1, 4, figsize=(7.1, 2.05), constrained_layout=True)
        for ax, family, title in zip(axes[:2], ['two_circles', 'clifford_torus'], ['(a) Two circles', '(b) Clifford torus']):
            r = [r for r in result['dictionary']['rows'] if r['family'] == family]
            ax.plot([x['atoms'] for x in r], [x['cover_radius'] / x['eta'] for x in r], '-o', ms=3, color=palette[0])
            good = [x for x in r if x['correct']]
            ax.scatter([x['atoms'] for x in good], [x['cover_radius'] / x['eta'] for x in good], marker='s', s=26, facecolors='none', edgecolors=palette[2], label='Correct rank')
            ax.axhline(1, ls='--', color='0.4', lw=.8, label='Certificate')
            ax.set(xscale='log', yscale='log', xlabel='Dictionary atoms', ylabel=r'$\varepsilon/\eta$', title=title)
            if family == 'two_circles': ax.legend(loc='upper right', frameon=False)
        ax = axes[2]
        for j, family in enumerate(['split_components', 'ellipse_deformation']):
            r = [r for r in result['mismatch']['rows'] if r['family'] == family]
            label = 'Split' if j == 0 else 'Ellipse'
            ax.plot([x['atoms'] for x in r], [x['distortion'] for x in r], '-o', ms=3, color=palette[j], label=label)
            ax.axhline(r[0]['reference_delta'], color=palette[j], ls='--', lw=.8)
        ax.set(xscale='log', xlabel='Dictionary atoms', ylabel='Bottleneck distortion', title='(c) Mismatch limit')
        ax.legend(frameon=False)
        r = [x for x in result['mismatch']['code_rows'] if x['K'] == 12 and x['trial'] == 0]
        axes[3].plot([x['m'] for x in r], [x['payload_bits'] for x in r], '-o', ms=3, color=palette[0], label='Joint integer')
        axes[3].plot([x['m'] for x in r], [x['coordinate_wise_bits'] for x in r], '--s', ms=3, color=palette[1], label='Per endpoint')
        axes[3].set(xlabel='Endpoint cells m', ylabel='Fixed payload bits', title='(d) Local sync, K=12')
        axes[3].legend(frameon=False)
        save(fig, 'dictionary_mismatch')

    if 'reparam' in result:
        fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.05), constrained_layout=True)
        rows = result['reparam']['rows']
        for ax, sigma in zip(axes[:2], [1.0, 2.0]):
            for j, (key, label) in enumerate([('naive', 'Naive'), ('raw_error', 'Noise reference'), ('debiased_error', 'Debiased')]):
                ds = [64, 256, 1024, 4096]
                groups = [[x[key] for x in rows if x['D_a'] == d and x['sigma'] == sigma] for d in ds]
                means = np.asarray([np.mean(g) for g in groups]); sd = np.asarray([np.std(g, ddof=1) / np.sqrt(len(g)) for g in groups])
                ax.errorbar(ds, means, yerr=1.96 * sd, color=palette[j], marker='o', markersize=3, lw=1, label=label, capsize=2)
            for key, color in [('raw_bound', palette[1]), ('debiased_bound', palette[2])]:
                bound = [next(x[key] for x in rows if x['D_a'] == d and x['sigma'] == sigma) for d in ds]
                ax.plot(ds, bound, ':', color=color, lw=1)
            ax.set(xscale='log', yscale='log', xlabel=r'Ambient dimension $D_a$', ylabel='Bottleneck distance', title=f'({"a" if sigma == 1 else "b"}) Total noise sigma={sigma:g}')
        axes[0].legend(frameon=False, loc='lower left')
        r = result['reparam']['calibration']; levels = sorted(set(x['relative_variance_error'] for x in r))
        axes[2].plot(levels, [np.mean([x['error'] for x in r if x['relative_variance_error'] == v]) for v in levels], '-o', ms=3, color=palette[2], label='Measured mean')
        axes[2].plot(levels, [next(x['bound'] for x in r if x['relative_variance_error'] == v) for v in levels], ':', color=palette[1], label='Calibrated bound')
        axes[2].set(xlabel='Relative variance error', ylabel='Bottleneck distance', title='(c) Calibration, D=1024')
        axes[2].legend(frameon=False)
        save(fig, 'noise_reparameterization')

    if 'channel' in result:
        fig, axes = plt.subplots(1, 4, figsize=(7.1, 2.1), constrained_layout=True)
        rows = result['channel']['probabilities']
        for j, Da in enumerate([64, 256, 1024, 4096]):
            for ax, event in zip(axes[:2], ['raw_spurious', 'partition_failure']):
                r = [x for x in rows if x['D_a'] == Da and x['event'] == event]
                ax.plot([x['sigma'] if event == 'raw_spurious' else x['normalized_sigma'] for x in r], [x['probability'] for x in r], '-', color=palette[j], label=str(Da), lw=1)
        axes[0].set(xlim=(.15, .8), xlabel='Total noise sigma', ylabel='Probability', title='(a) Raw fragmentation')
        axes[0].axvline(np.sqrt(2) / 4, color='0.5', ls=':', lw=.8)
        axes[0].legend(title=r'$D_a$', frameon=False, fontsize=6, title_fontsize=7)
        axes[1].set(xlim=(.15, .65), xlabel=r'$\sigma/[\Delta(D_a/\ln m)^{1/4}]$', ylabel='Probability', title='(b) Partition failure')
        ax = axes[2]
        r = [x for x in rows if x['D_a'] == 1024 and x['event'] == 'zero_cross_contact']
        r = [x for x in r if .2 <= x['normalized_sigma'] <= .65]
        xs = [x['normalized_sigma'] for x in r]
        ax.fill_between(xs, [x['lower_bound'] for x in r], [x['upper_bound'] for x in r], color=palette[0], alpha=.15, label='Proved bounds')
        ax.errorbar(xs, [x['probability'] for x in r], yerr=np.array([[x['probability'] - x['ci_low'] for x in r], [x['ci_high'] - x['probability'] for x in r]]), fmt='o-', ms=3, capsize=2, lw=1, color=palette[0], label='95% Wilson')
        ax.axvline(2 ** -1.25, color=palette[1], ls=':', lw=1)
        ax.set(xlabel='Normalized noise', ylabel='Probability', title='(c) Zero contact, D=1024')
        ax.legend(frameon=False, fontsize=6)
        ax = axes[3]; r = result['channel']['block_reliability']
        bracketed = [x for x in r if x['graph_last_below_sigma'] is not None
                     and x['graph_first_at_or_above_sigma'] is not None]
        lower = np.array([x['graph_last_below_sigma'] for x in bracketed])
        upper = np.array([x['graph_first_at_or_above_sigma'] for x in bracketed])
        center = np.sqrt(lower * upper)
        ax.errorbar([x['D_a'] for x in bracketed], center,
            yerr=np.array([center - lower, upper - center]), fmt='none',
            capsize=3, elinewidth=1.5, color=palette[0], label='Graph: grid bracket')
        ax.plot([x['D_a'] for x in r], [x['known_center_block_sigma'] for x in r],
            '--s', ms=3, color=palette[2], label='Known centers: exact')
        ax.set(xscale='log', yscale='log', xlabel=r'Ambient dimension $D_a$',
            ylabel='Total noise sigma', title='(d) 10% block crossing')
        ax.legend(frameon=False, fontsize=6)
        save(fig, 'channel_events')

    if 'task' in result and 'digits' in result:
        fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.05), constrained_layout=True)
        r = result['task']['quantization']
        axes[0].plot([x['step'] for x in r], [x['excess'] for x in r], '-o', ms=3, color=palette[0], label='Measured excess')
        axes[0].axhline(result['task']['exact_branch_entropy'], ls='--', color=palette[1], label='Exact limit')
        axes[0].set(xscale='log', xlabel='Score quantization step', ylabel='Extra entropy (bits)', title='(a) Branch-aware torus task')
        axes[0].legend(frameon=False)
        for j, pair in enumerate(['pixels64->pca16', 'pca16->pixels64']):
            r = result['digits']['dictionary'][pair]
            axes[1].plot(r['K'], r['d_topo'], '-o', ms=3, color=palette[j], label=pair.replace('pixels64', '64D').replace('pca16', '16D'))
        axes[1].axhline(r['floor'], color='0.5', ls=':', label='Full-cloud mismatch')
        axes[1].set(xscale='log', xlabel='Dictionary atoms', ylabel='Bottleneck distortion', title='(b) Held-out digit descriptors')
        axes[1].legend(frameon=False, fontsize=6)
        for j, model in enumerate(['pixels64', 'pca16']):
            r = result['digits']['cliff'][model]['rows']
            axes[2].plot([x['sigma'] for x in r], [x['kb_aware_acc'] for x in r], '-o', ms=3, color=palette[j], label='64 raw pixels' if j == 0 else '16 PCA coordinates')
        axes[2].set(xlabel='Total noise sigma', ylabel='Held-out class accuracy', ylim=(0,1), title='(c) Disjoint-reference receiver')
        axes[2].legend(frameon=False)
        save(fig, 'task_digits')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'code/results/topology')
    parser.add_argument('--stages', default='dictionary,reparam,channel,mismatch,task,digits')
    parser.add_argument('--reparam-trials', type=int, default=20)
    parser.add_argument('--channel-trials', type=int, default=120)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    stages = args.stages.split(',')
    if set(stages) - {'dictionary', 'reparam', 'channel', 'mismatch', 'task', 'digits'}:
        parser.error('unknown stage')
    start = time.perf_counter()
    report = {}
    functions = {'dictionary': dictionary, 'reparam': lambda out: reparameterization(out, args.reparam_trials),
        'channel': lambda out: channel(out, args.channel_trials), 'mismatch': mismatch_sync,
        'task': task, 'digits': digits}
    with threadpool_limits(limits=1):
        for name in stages:
            t0 = time.perf_counter()
            report[name] = functions[name](args.output)
            write_json(args.output / (name + '_summary.json'), report[name])
            print(f'{name}: completed in {time.perf_counter() - t0:.1f}s', flush=True)
        # Existing measurements from independent stage calls are eligible for
        # figure construction, but only executed stages enter this run report.
        available = dict(report)
        for name in functions:
            path = args.output / (name + '_summary.json')
            if name not in available and path.exists():
                available[name] = json.loads(path.read_text())
        figures(args.output, available)
    sources = [HERE, ROOT / 'code/topology/tda.py', ROOT / 'code/topology/reeb.py',
        ROOT / 'code/topology/sync.py', ROOT / 'code/embeddings/topology_exp.py',
        ROOT / 'code/embeddings/run_topology.py']
    metadata = {'seed': SEED, 'stages': stages, 'reparam_trials': args.reparam_trials,
        'channel_trials': args.channel_trials, 'python': platform.python_version(),
        'dependencies': {name: importlib.metadata.version(name) for name in
            ['numpy', 'scipy', 'ripser', 'persim', 'gudhi', 'scikit-learn', 'matplotlib']},
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'elapsed_s': time.perf_counter() - start,
        'monte_carlo_scope': 'probabilities are per fixed cloud/condition; no family-wise coverage is claimed',
        'radius_convention': 'maximum edge length divided by two',
        'noise_convention': 'per-axis variance sigma^2 / D_a',
        'figure_error_bars': 'reparam: 1.96 Monte Carlo standard errors; event probabilities: Wilson 95%'}
    write_json(args.output / 'results.json', {'metadata': metadata, 'results': report})
    print(f'Results and four vector figures: {args.output}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
