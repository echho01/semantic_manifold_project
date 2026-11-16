"""Analytic/synthetic checks for the proposed semantic-manifold theory.
No learned model, benchmark performance, or empirical communication gain is claimed.
Usage: python validate_theory.py --output .
Dependencies: numpy, scipy, matplotlib (Python >=3.10).
"""
from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rips_diagrams(points: np.ndarray) -> dict[int, np.ndarray]:
    """Exact Z2 boundary reduction of the complete 2-skeleton; H0/H1 only."""
    n = len(points)
    distance = np.linalg.norm(points[:, None] - points[None, :], axis=-1)
    simplices = [(0., (i,)) for i in range(n)]
    for size in (2, 3):
        for simplex in itertools.combinations(range(n), size):
            filtration = max(distance[a,b] for a,b in itertools.combinations(simplex,2))
            simplices.append((float(filtration), simplex))
    simplices.sort(key=lambda item: (item[0], len(item[1]), item[1]))
    idx = {s:j for j,(_,s) in enumerate(simplices)}
    pivot_columns: dict[int, set[int]] = {}
    diagrams = {0:[], 1:[]}
    for j, (death,simplex) in enumerate(simplices):
        column = set() if len(simplex)==1 else {idx[f] for f in itertools.combinations(simplex,len(simplex)-1)}
        while column and max(column) in pivot_columns:
            column ^= pivot_columns[max(column)]
        if column:
            p = max(column)
            pivot_columns[p] = column
            birth, face = simplices[p]
            dim = len(face)-1
            if dim in diagrams and death > birth + 1e-12:
                diagrams[dim].append((birth,death))
    return {dim:np.asarray(bars,dtype=float).reshape(-1,2) for dim,bars in diagrams.items()}


def bottleneck(a:np.ndarray,b:np.ndarray)->float:
    """Bottleneck metric with L-infinity cost and unlimited diagonal matches."""
    n,m=len(a),len(b)
    if n+m == 0: return 0.
    costs=np.full((n+m,m+n),np.inf)
    if n and m: costs[:n,:m]=np.max(np.abs(a[:,None,:]-b[None,:,:]),axis=2)
    for i in range(n): costs[i,m+i]=(a[i,1]-a[i,0])/2
    for j in range(m): costs[n+j,j]=(b[j,1]-b[j,0])/2
    costs[n:,m:]=0.
    candidates=np.unique(costs[np.isfinite(costs)])
    lo,hi=0,len(candidates)-1
    while lo<hi:
        mid=(lo+hi)//2
        adjacency=csr_matrix((costs <= candidates[mid]+1e-13).astype(np.int8))
        matching=maximum_bipartite_matching(adjacency,perm_type='column')
        if (matching>=0).all(): hi=mid
        else: lo=mid+1
    return float(candidates[lo])


def main(output:Path)->None:
    output.mkdir(parents=True,exist_ok=True)
    report={'status':'analytic_and_synthetic_checks_only','checks':{}}
    # Unit two-sphere: the exact homogeneous Gibbs RDF.
    rows=[]
    for s in [4,10,30,100,300,1000,3000]:
        z=2*math.pi*quad(lambda t:math.exp(-s*t*t)*math.sin(t),0,math.pi,epsabs=1e-14)[0]
        d=2*math.pi*quad(lambda t:t*t*math.exp(-s*t*t)*math.sin(t),0,math.pi,epsabs=1e-14)[0]/z
        r=math.log(4*math.pi)-math.log(z)-s*d
        base=math.log(4*math.pi)-math.log(math.pi*math.e*d)
        rows.append([s,d,r,base,(r-base)/d])
    sphere=np.array(rows)
    np.savetxt(output/'sphere_rd.csv',sphere,delimiter=',',header='s,D,R_nats,leading_R_nats,correction_over_D',comments='')
    assert abs(sphere[-1,-1]-1/6)<2e-5
    report['checks']['sphere']={'target':1/6,'measured_at_s_3000':float(sphere[-1,-1]),'pass':True}
    fig,ax=plt.subplots(figsize=(6.2,3.6))
    ax.semilogx(sphere[:,1],sphere[:,-1],'o-',label='Exact sphere RDF: (R - B) / D')
    ax.axhline(1/6,linestyle='--',label='Predicted coefficient: 1/6')
    ax.set_xlabel('Squared geodesic distortion D'); ax.set_ylabel('Finite-distortion coefficient')
    ax.legend(fontsize=9); ax.grid(alpha=.25); fig.tight_layout()
    fig.savefig(output/'sphere_check.png',dpi=220); fig.savefig(output/'sphere_check.pdf'); plt.close(fig)
    # Nonuniform density on S1, using the explicit forward Gibbs kernel.
    a=.5
    entropy=lambda amp: -quad(lambda t: ((1+amp*math.cos(t))/(2*math.pi))*math.log((1+amp*math.cos(t))/(2*math.pi)),-math.pi,math.pi,epsabs=1e-12)[0]
    h=entropy(a); fisher=1-math.sqrt(1-a*a)
    rows=[]
    for s in [10,30,100,300,1000,3000]:
        z=quad(lambda t:math.exp(-s*t*t),-math.pi,math.pi,epsabs=1e-14)[0]
        d=quad(lambda t:t*t*math.exp(-s*t*t),-math.pi,math.pi,epsabs=1e-14)[0]/z
        c=quad(lambda t:math.cos(t)*math.exp(-s*t*t),-math.pi,math.pi,epsabs=1e-14)[0]/z
        upper=entropy(a*c)-math.log(z)-s*d
        base=h-.5*math.log(2*math.pi*math.e*d)
        rows.append([s,d,upper,base,(upper-base)/d])
    circle=np.array(rows); np.savetxt(output/'nonuniform_circle.csv',circle,delimiter=',',header='s,D,achievable_I,leading_B,gap_over_D',comments='')
    assert abs(circle[-1,-1]-fisher/2)<2e-5
    report['checks']['nonuniform_circle']={'target_J_over_2':fisher/2,'measured':float(circle[-1,-1]),'pass':True}
    # Topology: same indexed vertices, arbitrary bounded perturbations.
    rng=np.random.default_rng(818)
    angles=np.linspace(0,2*math.pi,20,endpoint=False)
    points=np.c_[np.cos(angles),np.sin(angles)]
    source=rips_diagrams(points)
    topology=[]
    for scale in [.005,.01,.02,.04,.08,.15,.3]:
        for rep in range(4):
            noise=rng.normal(size=points.shape)*scale
            target=rips_diagrams(points+noise)
            eps=float(np.linalg.norm(noise,axis=1).max())
            dt=max(bottleneck(source[k],target[k]) for k in [0,1])
            topology.append([scale,rep,eps,dt,2*eps,dt/(2*eps)])
            assert dt<=2*eps+1e-9
    topology=np.array(topology); np.savetxt(output/'topology_stability.csv',topology,delimiter=',',header='scale,repeat,epsilon,bottleneck,certificate,ratio',comments='')
    report['checks']['topology']={'trials':len(topology),'max_ratio_to_bound':float(topology[:,-1].max()),'pass':True}
    fig,ax=plt.subplots(figsize=(6.2,3.6)); ax.scatter(topology[:,2],topology[:,3],label='Computed persistence distance')
    xx=np.array([0,topology[:,2].max()*1.05]); ax.plot(xx,2*xx,'--',label='Certificate: 2 epsilon')
    ax.set_xlabel('Maximum vertex displacement epsilon'); ax.set_ylabel('Bottleneck distance')
    ax.legend(fontsize=9); ax.grid(alpha=.25);fig.tight_layout();fig.savefig(output/'topology_check.png',dpi=220);fig.savefig(output/'topology_check.pdf');plt.close(fig)
    # Exact minimizing movement for a discrete quadratic bending surrogate.
    n=64; grid=np.linspace(0,2*np.pi,n,endpoint=False)
    f0=np.sin(grid)+.3*np.sin(10*grid); f=f0.copy()
    L=2*np.eye(n)-np.roll(np.eye(n),1,0)-np.roll(np.eye(n),-1,0)
    la,lb,lf,tau=.5,4.,.2,.4
    A=la*L+lb*(L@L)+lf*np.eye(n)
    energy=lambda v: float(.5*la*v@L@v+.5*lb*np.linalg.norm(L@v)**2+.5*lf*np.linalg.norm(v-f0)**2)
    energies=[energy(f)]
    for _ in range(40):
        f=np.linalg.solve(np.eye(n)+tau*A,f+tau*lf*f0); energies.append(energy(f))
    assert max(np.diff(energies))<1e-10
    np.savetxt(output/'variational_energy.csv',np.c_[np.arange(len(energies)),energies],delimiter=',',header='iteration,energy',comments='')
    report['checks']['variational_step']={'initial_energy':energies[0],'final_energy':energies[-1],'monotone':True,'pass':True}
    (output/'results.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    (output/'REPORT.md').write_text('# 解析与合成算例核验\n\n本目录核验推导中的常数、渐近阶和离散稳定性；数据由给定数学模型生成。\n\n'+ '\n'.join(f'## {name}\n```json\n{json.dumps(value,ensure_ascii=False,indent=2)}\n```\n' for name,value in report['checks'].items())+'\n真实图像/点云语义模型、信道实测及与TopoJSCC的训练对照列入后续实验。\n',encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path(__file__).parent)
    main(parser.parse_args().output)
