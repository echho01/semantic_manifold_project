#!/usr/bin/env python3
"""Reproducible mathematical checks; NOT learned-codec performance experiments.
Run: python scripts/validate_implementation.py --out results/implementation
Dependencies: numpy, scipy, matplotlib. Seed and numerical errors are recorded.
"""
import argparse, itertools, json, platform
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import linear_sum_assignment, minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sphere(s):
    upper = min(40., np.pi*np.sqrt(s))
    f = lambda t: np.exp(-t*t)*np.sin(t/np.sqrt(s))/np.sqrt(s)
    z = 2*np.pi*quad(f, 0, upper, epsabs=1e-14, epsrel=1e-12)[0]
    d = 2*np.pi*quad(lambda t: t*t/s*f(t), 0, upper, epsabs=1e-14, epsrel=1e-12)[0]/z
    rate = np.log(4*np.pi)-np.log(z)-s*d
    base = np.log(4*np.pi)-np.log(np.pi*np.e*d)
    return dict(s=s, distortion=d, rate_nats=rate, correction=(rate-base)/d)


def rips(points):
    """Finite H0,H1 intervals of the COMPLETE Rips two-skeleton over F2."""
    n = len(points)
    distances = np.linalg.norm(points[:, None]-points[None, :], axis=-1)
    cells = [(0., (i,)) for i in range(n)]
    for size in (2,3):
        for simplex in itertools.combinations(range(n), size):
            cells.append((max(distances[i,j] for i,j in itertools.combinations(simplex,2)), simplex))
    cells.sort(key=lambda a: (a[0],len(a[1]),a[1]))
    index = {simplex:i for i,(_,simplex) in enumerate(cells)}
    pivots = {}; diagrams = [[],[]]
    for j,(death,simplex) in enumerate(cells):
        column = set() if len(simplex)==1 else {index[f] for f in itertools.combinations(simplex,len(simplex)-1)}
        while column and max(column) in pivots:
            column ^= pivots[max(column)]
        if column:
            birth_index=max(column); pivots[birth_index]=column
            birth,birth_cell=cells[birth_index]; dim=len(birth_cell)-1
            if dim<=1 and death>birth+1e-12:
                diagrams[dim].append((birth,death))
    return [np.array(d).reshape(-1,2) for d in diagrams]


def bottleneck(a,b):
    n,m=len(a),len(b)
    if n+m==0: return 0.
    cost=np.full((n+m,n+m),np.inf)
    if n and m: cost[:n,:m]=np.max(np.abs(a[:,None]-b[None,:]),axis=-1)
    for i in range(n):cost[i,m+i]=(a[i,1]-a[i,0])/2
    for j in range(m):cost[n+j,j]=(b[j,1]-b[j,0])/2
    cost[n:,m:]=0
    candidates=np.unique(cost[np.isfinite(cost)])
    lo,hi=0,len(candidates)-1
    while lo<hi:
        mid=(lo+hi)//2
        row,col=linear_sum_assignment((cost>candidates[mid]).astype(float))
        if np.all(cost[row,col]<=candidates[mid]):hi=mid
        else:lo=mid+1
    return float(candidates[lo])


def mutual(q,w):
    marginal=w@q
    return float(np.sum(w[:,None]*q*np.log(q/marginal[None,:])))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='results/implementation');args=parser.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20260915)
    report={'kind':'executed mathematical validation; no neural training','seed':20260915,'python':platform.python_version()}
    report['sphere']=[sphere(s) for s in [30.,100.,300.,1000.,3000.]]
    n=8192; angle=2*np.pi*np.arange(n)/n; dx=2*np.pi/n
    distance=np.minimum(angle,2*np.pi-angle);p=(1+.5*np.cos(angle))/(2*np.pi)
    h0=-np.sum(p*np.log(p))*dx;circle=[]
    for s in [30.,100.,300.,1000.]:
        k=np.exp(-s*distance**2);z=np.sum(k)*dx;k/=z
        pout=np.fft.irfft(np.fft.rfft(p)*np.fft.rfft(k),n=n)*dx
        d=np.sum(k*distance**2)*dx
        info=-np.sum(pout*np.log(pout))*dx-s*d-np.log(z)
        base=h0-.5*np.log(2*np.pi*np.e*d)
        circle.append({'s':s,'distortion':d,'correction':(info-base)/d})
    report['circle']={'expected_forward_coefficient':(1-np.sqrt(.75))/2,'observed':circle}
    points=rng.normal(size=(9,3));velocity=rng.normal(size=points.shape);v=np.exp(rng.normal(size=9));t=rng.normal(size=9)
    w=rng.uniform(.2,1,size=9);w/=w.sum();s=.7;lam=.3
    cg=np.sum((points[:,None]-points[None,:])**2,axis=-1)
    def kernel(theta):
        x=points+theta*velocity;cost=np.sum((x[:,None]-x[None,:])**2,axis=-1)
        mass=np.exp(-s*cost)*v[None,:]*np.exp(theta*t)[None,:]
        return mass/mass.sum(axis=1,keepdims=True)
    q=kernel(0);dc=2*np.sum((points[:,None]-points[None,:])*(velocity[:,None]-velocity[None,:]),axis=-1)
    a=t[None,:]-s*dc;f=np.log(q/(w@q)[None,:])+lam*cg
    grad=np.sum(w*(np.sum(q*f*a,axis=1)-np.sum(q*f,axis=1)*np.sum(q*a,axis=1)))
    fun=lambda theta:mutual(kernel(theta),w)+lam*np.sum(w[:,None]*kernel(theta)*cg)
    step=1e-6;fd=(fun(step)-fun(-step))/(2*step)
    scaled=np.exp(-(s/9)*(9*cg))*(9*v)[None,:];scaled/=scaled.sum(axis=1,keepdims=True)
    report['kernel']={'gradient':float(grad),'finite_difference':float(fd),'absolute_error':float(abs(grad-fd)),
                      'dilation_max_transition_error':float(np.max(abs(q-scaled)))}
    theta=np.linspace(0,2*np.pi,12,endpoint=False);clean=np.c_[np.cos(theta),np.sin(theta)]
    cleanD=rips(clean);cases=[]
    for sigma in [.005,.02,.05,.1]:
        for _ in range(10):
            noisy=clean+sigma*rng.normal(size=clean.shape);d=rips(noisy)
            eps=np.max(abs(np.linalg.norm(clean[:,None]-clean[None,:],axis=-1)-np.linalg.norm(noisy[:,None]-noisy[None,:],axis=-1)))
            db=max(bottleneck(cleanD[k],d[k]) for k in (0,1))
            assert db<=eps+1e-9
            cases.append({'sigma':sigma,'metric_bound':float(eps),'bottleneck':db})
    translated=rips(clean+np.array([7.,-4.]));shift=max(bottleneck(cleanD[k],translated[k]) for k in (0,1))
    report['persistence']={'trials':cases,'max_bound_ratio':max(x['bottleneck']/x['metric_bound'] for x in cases),'translation_distance':shift}
    c=np.array([.2,1.,3.,7.]);T=.8;power=np.sqrt(c)*np.sqrt(c).sum()/T
    solved=minimize(lambda p:p.sum(),power*1.5,method='SLSQP',bounds=[(1e-6,None)]*4,
                    constraints={'type':'ineq','fun':lambda p:T-np.sum(c/p)},options={'ftol':1e-11,'maxiter':1000})
    report['power']={'closed_form':power.tolist(),'constraint':float(np.sum(c/power)),
                     'optimizer_success':bool(solved.success),'relative_objective_error':float(abs(solved.fun-power.sum())/power.sum())}
    N,dim,alpha,total=12,4,.05,20000;variances=np.array([.003,.008,.02,.06]);tt=np.log(N*(N-1)/(2*alpha))
    radius=np.sqrt(2*variances.sum())+2*np.sqrt(variances.max()*tt)
    pairs=np.array(list(itertools.combinations(range(N),2)));fail=0
    for start in range(0,total,500):
        e=rng.normal(size=(min(500,total-start),N,dim))*np.sqrt(variances)
        maxerr=np.linalg.norm(e[:,pairs[:,0]]-e[:,pairs[:,1]],axis=-1).max(axis=1)
        fail+=int((maxerr>radius).sum())
    from scipy.stats import beta
    upper=1. if fail==total else float(beta.ppf(.95,fail+1,total-fail))
    report['gaussian']={'draws':total,'nominal_alpha':alpha,'radius':float(radius),'violations':fail,'one_sided_95pct_binomial_upper':upper}
    # Binary type B and binary shared map C. Both code endpoints know C.
    joint=np.array([[.4,.1],[.1,.4]]);pc=joint.sum(axis=0);cond=joint/pc[None,:];model=np.array([[.65,.25],[.35,.75]])
    entropy=-np.sum(joint*np.log(cond));cross=-np.sum(joint*np.log(model));kl=np.sum(joint*np.log(cond/model))
    report['knowledge']={'H_B_given_C_nats':float(entropy),'mismatch_KL_nats':float(kl),
                         'cross_entropy_nats':float(cross),'identity_error':float(abs(cross-entropy-kl))}
    assert report['kernel']['absolute_error']<1e-7
    assert report['kernel']['dilation_max_transition_error']<1e-12
    assert report['power']['relative_objective_error']<1e-6
    assert shift<1e-10
    (out/'validation.json').write_text(json.dumps(report,indent=2))
    # One chart per figure; standard Matplotlib styling.
    fig,ax=plt.subplots(figsize=(5.4,3.6))
    ax.plot([x['distortion'] for x in report['sphere']],[x['correction'] for x in report['sphere']],'o-',label='Sphere quadrature')
    ax.axhline(1/6,linestyle='--',label='Predicted 1/6');ax.set_xscale('log');ax.set_xlabel('Distortion D');ax.set_ylabel('(R - leading term) / D');ax.legend();fig.tight_layout();fig.savefig(out/'geometry_validation.pdf');fig.savefig(out/'geometry_validation.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.4,3.6))
    ax.scatter([x['metric_bound'] for x in cases],[x['bottleneck'] for x in cases],s=22,label='Executed perturbations')
    hi=max(x['metric_bound'] for x in cases)*1.04;ax.plot([0,hi],[0,hi],'--',label='Proved upper bound');ax.set_xlabel('Maximum pair-distance change');ax.set_ylabel('Bottleneck distance');ax.legend();fig.tight_layout();fig.savefig(out/'topology_validation.pdf');fig.savefig(out/'topology_validation.png',dpi=180);plt.close(fig)
    summary={k:v for k,v in report.items() if k not in ('persistence','sphere','circle')};summary['sphere_limit_check']=report['sphere'][-1];summary['persistence_max_ratio']=report['persistence']['max_bound_ratio']
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
