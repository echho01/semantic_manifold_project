#!/usr/bin/env python3
"""Fixed-cloud persistence, quantization and AWGN experiments for Paper 2.

The Rips parameter is the edge length. All channel trials perturb the same
clean point cloud. Finite H0 bars and essential H0 are recorded separately.
Run: python code/synthetic/run_paper2.py --output code/results/paper2
"""
from __future__ import annotations
import argparse
import csv
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import sys
import time
os.environ.setdefault("MPLCONFIGDIR", "/tmp/semantic-manifold-matplotlib")
import numpy as np
from scipy.stats import norm
from ripser import ripser
from persim import bottleneck
from threadpoolctl import threadpool_limits


def clean_circles(n_per_circle=64, separation=.6):
    theta = np.arange(n_per_circle) * 2 * np.pi / n_per_circle
    circle = np.column_stack([np.cos(theta), np.sin(theta)])
    center = 1 + separation / 2
    return np.concatenate([circle - [center, 0], circle + [center, 0]])


def finite_diagrams(points):
    diagrams = ripser(points, maxdim=1)["dgms"]
    essential = [int(np.sum(~np.isfinite(d[:,1]))) for d in diagrams]
    return [d[np.isfinite(d[:,1])] for d in diagrams], essential


def margin(diagram, theta):
    life = diagram[:,1] - diagram[:,0]
    return float(min(theta, np.min(np.abs(life-theta), initial=np.inf)))


def thresholds_for_signals(diagrams, counts=(1,2)):
    thresholds = []
    for diagram, count in zip(diagrams, counts):
        life = np.sort(diagram[:,1]-diagram[:,0])[::-1]
        v = float(life[count-1])
        u = float(life[count]) if len(life)>count else 0.
        if v <= u:
            raise ValueError("requested signal bars have no persistence gap")
        thresholds.append((u+v)/2)
    return thresholds


def counts(diagrams, theta):
    return [int(np.count_nonzero(d[:,1]-d[:,0]>th)) for d,th in zip(diagrams,theta)]


def quantize(points, bits, half_range):
    """Uniform midrise quantizer on [-L,L], b bits per coordinate."""
    if np.any(np.abs(points)>half_range):
        raise ValueError("source coordinates exceed the specified quantizer range")
    levels = 2**bits
    index = np.minimum(np.floor((points+half_range)/(2*half_range)*levels),levels-1)
    return -half_range+(index+.5)*(2*half_range/levels)


def wilson(failures, n):
    p = failures/n; z = norm.ppf(.975); den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z/den*np.sqrt(p*(1-p)/n+z*z/(4*n*n))
    return max(0.,center-half),min(1.,center+half)


def minimum_bits(point_budget, sigma, n, ambient, half_range, delta):
    noise_radius = sigma*(np.sqrt(ambient)+np.sqrt(2*np.log(n/delta)))
    remainder = point_budget-noise_radius
    if remainder<=0:
        return None
    # Strict inequality q_b+a_delta < point_budget.
    return max(1,int(np.floor(np.log2(np.sqrt(ambient)*half_range/remainder)))+1)


def run(output, trials=40, seed=20261019):
    started=time.time(); output.mkdir(parents=True,exist_ok=True)
    X=clean_circles(); n,ambient=X.shape; half_range=2.5; delta=.05
    clean,essential=finite_diagrams(X)
    schemes={"fixed": [.2,.2], "midpoint":thresholds_for_signals(clean)}
    budgets={key:min(margin(d,theta) for d,theta in zip(clean,value))/4
             for key,value in schemes.items()}
    expected={key:counts(clean,value) for key,value in schemes.items()}
    result={"metadata":{"schema_version":1,"seed":seed,"trials_per_condition":trials,
            "n_points":n,"ambient_dimension":ambient,"rips_parameter":"edge_length",
            "circle_radius":1.,"circle_separation":.6,"quantizer_half_range":half_range,
            "delta":delta,"python":sys.version.split()[0],"numpy":np.__version__,
            "ripser":version("ripser"),"persim":version("persim"),"platform":platform.platform(),
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "command":" ".join(sys.argv)},"clean_diagrams":[d.tolist() for d in clean],
            "essential_counts":essential,"thresholds":schemes,
            "expected_finite_counts":expected,"point_error_budgets":budgets,
            "conditions":[],"trials":[],"bit_certificates":[]}
    noise_sigmas=(0.,.002,.004,.008,.012,.02,.04,.08,.12,.20)
    conditions=[("noise",s,None) for s in noise_sigmas]+[("bits",.002,b) for b in range(2,11)]
    for condition_index,(kind,sigma,bits) in enumerate(conditions):
        quantized=X if bits is None else quantize(X,bits,half_range)
        qbound=0. if bits is None else np.sqrt(ambient)*half_range*2.**(-bits)
        abound=sigma*(np.sqrt(ambient)+np.sqrt(2*np.log(n/delta)))
        entries=[]
        for trial in range(trials):
            # Common noise for both threshold policies in every trial.
            # Each condition has its own deterministically identified stream.
            stream=np.random.SeedSequence([seed,condition_index,trial])
            rng=np.random.default_rng(stream)
            noise=rng.normal(scale=sigma,size=X.shape)
            Y=quantized+noise
            dgms,ess=finite_diagrams(Y)
            max_error=float(np.linalg.norm(Y-X,axis=1).max())
            distances=[float(bottleneck(a,b)) for a,b in zip(clean,dgms)]
            if max(distances)>2*max_error+2e-6:
                raise AssertionError("Rips stability check failed")
            entry={"condition":kind,"sigma_axis":sigma,"bits_per_coordinate":bits,
                   "trial":trial,"max_point_error":max_error,
                   "max_noise_norm":float(np.linalg.norm(noise,axis=1).max()),
                   "bottleneck_h0":distances[0],"bottleneck_h1":distances[1]}
            for scheme,theta in schemes.items():
                observed=counts(dgms,theta)
                failure=observed!=expected[scheme] or ess!=essential
                certified=max_error<budgets[scheme]
                if certified and failure:
                    raise AssertionError("count certificate violated")
                entry[scheme+"_finite_counts"]=observed
                entry[scheme+"_outage"]=int(failure)
                entry[scheme+"_sample_certified"]=bool(certified)
            entries.append(entry)
        row={"condition":kind,"sigma_axis":sigma,"bits_per_coordinate":bits,
             "trials":trials,"q_bound":float(qbound),"gaussian_radius_bound":float(abound),
             "mean_bottleneck_h0":float(np.mean([e["bottleneck_h0"] for e in entries])),
             "mean_bottleneck_h1":float(np.mean([e["bottleneck_h1"] for e in entries])),
             "mean_max_point_error":float(np.mean([e["max_point_error"] for e in entries]))}
        for scheme in schemes:
            failures=sum(e[scheme+"_outage"] for e in entries)
            lo,hi=wilson(failures,trials)
            row.update({scheme+"_outages":failures,scheme+"_outage_rate":failures/trials,
                        scheme+"_ci95_low":lo,scheme+"_ci95_high":hi,
                        scheme+"_probability_certified":bool(qbound+abound<budgets[scheme]),
                        scheme+"_sample_certified_fraction":float(np.mean([e[scheme+"_sample_certified"] for e in entries]))})
        result["conditions"].append(row);result["trials"].extend(entries)
        print(f"{kind}: sigma={sigma:.3f}, b={bits}, outage fixed={row['fixed_outage_rate']:.3f}, midpoint={row['midpoint_outage_rate']:.3f}",flush=True)
    for sigma in (0.,.002,.004,.006,.008,.010,.012):
        row={"sigma_axis":sigma}
        for scheme,budget in budgets.items():
            b=minimum_bits(budget,sigma,n,ambient,half_range,delta)
            row[scheme+"_minimum_bits_per_coordinate"]=b
            row[scheme+"_payload_bits"]=None if b is None else n*ambient*b
        result["bit_certificates"].append(row)
    result["checks"]={"stability_violations":0,"count_certificate_violations":0,
        "total_trials":len(result["trials"]),
        "midpoint_budget_over_fixed":budgets["midpoint"]/budgets["fixed"]}
    result["metadata"]["elapsed_seconds"]=time.time()-started
    (output/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    np.savez_compressed(output/"clean_cloud.npz",points=X,diagram_h0=clean[0],diagram_h1=clean[1])
    for key in ("conditions","bit_certificates"):
        with (output/(key+".csv")).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(result[key][0]));writer.writeheader();writer.writerows(result[key])
    plot(result,output)
    print(json.dumps({"budgets":budgets,"checks":result["checks"],"bit_certificates":result["bit_certificates"]},indent=2))
    return result


def plot(result,output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"legend.fontsize":7,
        "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none","legend.frameon":False})
    fig,axes=plt.subplots(1,2,figsize=(7.16,2.7),layout="constrained")
    for ax,kind,xkey,title in zip(axes,("noise","bits"),("sigma_axis","bits_per_coordinate"),
                                 ("(a) AWGN with exact coordinates","(b) Quantized coordinates, noise std. 0.002")):
        rows=[r for r in result["conditions"] if r["condition"]==kind]
        x=[r[xkey] for r in rows]
        for scheme,color,marker,label in (("fixed","#D55E00","s","Fixed threshold 0.2"),
                                          ("midpoint","#0072B2","o","Midpoint thresholds")):
            y=[r[scheme+"_outage_rate"] for r in rows]
            ax.plot(x,y,color=color,marker=marker,markersize=4,label=label)
            ax.fill_between(x,[r[scheme+"_ci95_low"] for r in rows],[r[scheme+"_ci95_high"] for r in rows],color=color,alpha=.12)
        ax.set(title=title,ylabel="Feature-count outage probability",ylim=(-.02,1.02));ax.grid(alpha=.18)
    axes[0].set_xlabel("Per-axis noise standard deviation")
    axes[1].set_xlabel("Quantizer bits per coordinate")
    axes[0].legend(loc="lower right", bbox_to_anchor=(1, .10))
    for ext in ("pdf","svg","png"):
        fig.savefig(output/f"persistence_reliability.{ext}",dpi=220)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1]/"results"/"paper2")
    parser.add_argument("--trials",type=int,default=40)
    parser.add_argument("--seed",type=int,default=20261019)
    args=parser.parse_args()
    with threadpool_limits(limits=1):run(args.output,args.trials,args.seed)
