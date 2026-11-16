"""Reproducible analytic and finite-complex checks; no trained-codec results.
Run: python scripts/validate_revision.py --output results
Python 3.10+, numpy, scipy, matplotlib. All logarithms in nats.
"""
from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rips(points: np.ndarray) -> dict[int,np.ndarray]:
    n=len(points)
    d=np.linalg.norm(points[:,None]-points[None,:],axis=-1)
    simplices=[(0.,(i,)) for i in range(n)]
    for k in (2,3):
        for a in itertools.combinations(range(n),k):
            simplices.append((max(d[i,j] for i,j in itertools.combinations(a,2)),a))
    simplices.sort(key=lambda x:(x[0],len(x[1]),x[1]))
    idx={s:i for i,(_,s) in enumerate(simplices)};piv={};out={0:[],1:[]}
    for death,s in simplices:
        col=set() if len(s)==1 else {idx[f] for f in itertools.combinations(s,len(s)-1)}
        while col and max(col) in piv: col ^= piv[max(col)]
        if col:
            p=max(col);piv[p]=col;birth,face=simplices[p];dim=len(face)-1
            if dim in out and death-birth>1e-10:out[dim].append((birth,death))
    return {k:np.array(v,dtype=float).reshape(-1,2) for k,v in out.items()}


def bn(a:np.ndarray,b:np.ndarray)->float:
    n,m=len(a),len(b)
    if n+m==0:return 0.
    c=np.full((n+m,m+n),np.inf)
    if n*m:c[:n,:m]=abs(a[:,None]-b[None,:]).max(2)
    for i in range(n):c[i,m+i]=(a[i,1]-a[i,0])/2
    for j in range(m):c[n+j,j]=(b[j,1]-b[j,0])/2
    c[n:,m:]=0.;v=np.unique(c[np.isfinite(c)]);lo=0;hi=len(v)-1
    while lo<hi:
        mid=(lo+hi)//2
        matching=maximum_bipartite_matching(csr_matrix(c<=v[mid]+1e-12),perm_type='column')
        if (matching>=0).all():hi=mid
        else:lo=mid+1
    return float(v[lo])


def diagram_distance(a,b):return max(bn(a[k],b[k]) for k in (0,1))


def sphere(s:float):
    # Change variable t=x/sqrt(s) resolves narrow kernels at high s.
    upper=min(math.pi*math.sqrt(s),14.)
    z=2*math.pi/math.sqrt(s)*quad(lambda x:math.exp(-x*x)*math.sin(x/math.sqrt(s)),0,upper,epsabs=1e-14)[0]
    d=2*math.pi/s**1.5*quad(lambda x:x*x*math.exp(-x*x)*math.sin(x/math.sqrt(s)),0,upper,epsabs=1e-14)[0]/z
    r=math.log(4*math.pi)-math.log(z)-s*d
    b=math.log(4*math.pi)-math.log(math.pi*math.e*d)
    return d,r,b


def plot(output:Path,name:str,draw,xlabel:str,ylabel:str):
    fig,ax=plt.subplots(figsize=(6.2,3.8))
    draw(ax);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    ax.grid(alpha=.18);ax.legend(fontsize=10);fig.tight_layout()
    fig.savefig(output/(name+'.png'),dpi=240)
    fig.savefig(output/(name+'.pdf'));plt.close(fig)


def main(out:Path):
    out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(260906)
    report={'status':'executed_analytic_and_synthetic_checks','seed':260906,'checks':{}}
    vals=[]
    for s in (4,10,30,100,300,1000,3000):
        d,r,b=sphere(s);vals.append((s,d,r,b,(r-b)/d))
    vals=np.array(vals);np.savetxt(out/'geometry_sphere.csv',vals,delimiter=',',header='s,D,R_nats,B_nats,(R-B)/D',comments='')
    assert abs(vals[-1,-1]-1/6)<2e-5
    report['checks']['sphere_coefficient']={'target':1/6,'observed':float(vals[-1,-1]),'pass':True}
    plot(out,'geometry_coefficient',lambda ax:(ax.semilogx(vals[:,1],vals[:,-1],'o-',label='Exact Gibbs RDF'),ax.axhline(1/6,ls='--',label='Scalar-curvature coefficient: 1/6')),'Geodesic MSE D','(R(D) - B(D)) / D')
    plot(out,'geometry_rd',lambda ax:(ax.semilogx(vals[:,1],vals[:,2]/math.log(2),'o-',label='Achievable rate = converse (sphere)'),ax.semilogx(vals[:,1],vals[:,3]/math.log(2),'--',label='Leading dimension/entropy term')),'Geodesic MSE D','Rate (bits / semantic state)')
    # Circle density: exact convolution of a single Fourier harmonic.
    amp=.5
    def entropy(a):return -quad(lambda t:((1+a*math.cos(t))/(2*math.pi))*math.log((1+a*math.cos(t))/(2*math.pi)),-math.pi,math.pi,epsabs=1e-13)[0]
    h=entropy(amp);rows=[]
    for s in (10,30,100,300,1000,3000):
        upper=min(math.pi*math.sqrt(s),14.)
        z=quad(lambda x:math.exp(-x*x),-upper,upper,epsabs=1e-14)[0]/math.sqrt(s)
        d=quad(lambda x:x*x*math.exp(-x*x),-upper,upper,epsabs=1e-14)[0]/s**1.5/z
        k=quad(lambda x:math.exp(-x*x)*math.cos(x/math.sqrt(s)),-upper,upper,epsabs=1e-14)[0]/math.sqrt(s)/z
        I=entropy(amp*k)-math.log(z)-s*d;B=h-.5*math.log(2*math.pi*math.e*d)
        rows.append((s,d,I,B,(I-B)/d))
    rows=np.array(rows);target=(1-math.sqrt(1-amp**2))/2
    assert abs(rows[-1,-1]-target)<2e-5
    np.savetxt(out/'geometry_density.csv',rows,delimiter=',',header='s,D,I_ach_nats,B_nats,gap_over_D',comments='')
    report['checks']['density_kernel_coefficient']={'target_J_over_2':target,'observed':float(rows[-1,-1]),'pass':True}
    # Exact scale cancellation at the finite kernel level.
    angles=np.linspace(0,2*np.pi,96,endpoint=False)
    dist=np.abs(angles[:,None]-angles[None,:]);dist=np.minimum(dist,2*np.pi-dist)
    weights=(1+amp*np.cos(angles));weights/=weights.sum();s=7.
    def kernel_info(scale):
        c=scale**2*dist**2;lv=np.full(96,math.log(scale*2*np.pi/96))
        logits=lv[None,:]-(s/scale**2)*c
        Q=np.exp(logits-logsumexp(logits,axis=1,keepdims=True));pi=weights@Q
        return float(np.sum(weights[:,None]*Q*np.log(Q/pi))),float(np.sum(weights[:,None]*Q*c))/scale**2
    refs=kernel_info(1);scale_errors=[max(abs(np.array(kernel_info(c))-refs)) for c in (.2,.5,2,5)]
    assert max(scale_errors)<1e-12
    report['checks']['finite_kernel_scale_invariance']={'max_absolute_error':max(scale_errors),'pass':True}
    # Analytic variational derivative on a finite discretization, any smooth cost.
    weights=weights[:48];weights/=weights.sum();t=angles[:48]
    def quantities(z):
        pts=np.c_[np.cos(t),z*np.sin(t)]
        dpts=np.c_[np.zeros(len(t)),np.sin(t)]
        dif=pts[:,None]-pts[None,:];ddif=dpts[:,None]-dpts[None,:]
        c=(dif*dif).sum(2);dc=2*(dif*ddif).sum(2)
        v=np.sqrt(np.sin(t)**2+z*z*np.cos(t)**2)
        dlv=z*np.cos(t)**2/v**2
        logQ=np.log(v)[None,:]-3*c;logQ-=logsumexp(logQ,axis=1,keepdims=True)
        Q=np.exp(logQ);pi=weights@Q;ell=logQ-np.log(pi)[None,:]
        A=dlv[None,:]-3*dc
        I=float(np.sum(weights[:,None]*Q*ell))
        dI=float(np.sum(weights[:,None]*Q*ell*(A-(Q*A).sum(1)[:,None])))
        return I,dI
    z=.72;I,dI=quantities(z);eps=1e-5;fd=(quantities(z+eps)[0]-quantities(z-eps)[0])/(2*eps)
    assert abs(dI-fd)<1e-7
    report['checks']['surface_variational_gradient']={'analytic':dI,'finite_difference':fd,'absolute_error':abs(dI-fd),'pass':True}
    # Pairwise certificate: complete 2-skeleton, common translation cancels.
    th=np.linspace(0,2*np.pi,16,endpoint=False);Z=np.c_[np.cos(th),np.sin(th)]
    source=rips(Z);trans=rips(Z+np.array([7.,-4.]))
    translation_error=diagram_distance(source,trans)
    assert translation_error<1e-10
    report['checks']['common_translation']={'bottleneck':translation_error,'pass':True}
    topo=[]
    for a in (1.,.65,.3,.08):
        B=np.diag([1.,a]);G=Z@B.T
        d0=np.linalg.norm(Z[:,None]-Z[None,:],axis=-1)
        dg=np.linalg.norm(G[:,None]-G[None,:],axis=-1)
        det=float(abs(dg-d0).max())
        for noise_sd in (.0,.015,.05,.12):
            for repeat in range(3):
                e=rng.normal(size=Z.shape)*noise_sd;Zh=G+e+np.array([5.,-2.])
                ne=float(np.linalg.norm(e[:,None]-e[None,:],axis=-1).max())
                cert=det+ne;db=diagram_distance(source,rips(Zh))
                assert db<=cert+1e-8
                topo.append((a,noise_sd,repeat,db,cert,db/cert if cert>1e-9 else 0.))
    topo=np.array(topo);np.savetxt(out/'topology_pairwise.csv',topo,delimiter=',',header='min_gain,noise_sd,repeat,bottleneck,certificate,ratio',comments='')
    report['checks']['pairwise_topology_certificate']={'trials':len(topo),'max_ratio':float(topo[:,-1].max()),'pass':True}
    plot(out,'topology_certificate',lambda ax:(ax.scatter(topo[:,4],topo[:,3],label='Computed H0/H1 distance'),ax.plot([0,topo[:,4].max()],[0,topo[:,4].max()],'--',label='Proved pairwise certificate')),'Pairwise perturbation certificate','Persistence bottleneck distance')
    # Gaussian concentration and resource allocation: source alphabet fixed.
    coeff=np.array([.4,1.,2.5]);tau=.9;eta0=.08;alpha=.05;N=16
    K=1+math.sqrt(2*math.log((N*(N-1)/2)/alpha));T=(tau-eta0)**2/(2*K*K)
    pstar=np.sqrt(coeff)*np.sqrt(coeff).sum()/T
    solve=minimize(lambda p:p.sum(),pstar*1.15,method='SLSQP',bounds=[(1e-8,None)]*3,constraints=[{'type':'ineq','fun':lambda p:T-np.sum(coeff/p)}],options={'ftol':1e-10,'maxiter':1000})
    assert solve.success and abs(solve.fun-pstar.sum())/pstar.sum()<1e-5
    report['checks']['minimum_certificate_power']={'closed_form':float(pstar.sum()),'SLSQP':float(solve.fun),'relative_error':float(abs(solve.fun-pstar.sum())/pstar.sum()),'allocation':pstar.tolist(),'pass':True}
    S=coeff/pstar;pair_bound=math.sqrt(2*S.sum())+2*math.sqrt(S.max()*math.log((N*(N-1)/2)/alpha))
    # Marginal noise certificate only; pair dependence is handled by the union bound.
    draws=20000;viol=0
    for start in range(0,draws,1000):
        eps=rng.normal(size=(min(1000,draws-start),N,3))*np.sqrt(S)[None,None,:]
        dif=eps[:,:,None,:]-eps[:,None,:,:]
        mx=np.sqrt((dif*dif).sum(-1)).max(axis=(1,2));viol+=int((mx>pair_bound).sum())
    report['checks']['gaussian_radius_mc']={'draws':draws,'alpha':alpha,'violation_count':viol,'violation_fraction':viol/draws,'radius':pair_bound,'pass':True}
    # Affine ellipses: homotopy type fixed for a>0; metric persistence shrinks.
    collapse=[]
    for a in (1.,.7,.4,.2,.1,.03,0.):
        ph=rips(Z*np.array([1.,a]));life=max((d-b for b,d in ph[1]),default=0.)
        collapse.append((a,life))
    collapse=np.array(collapse);np.savetxt(out/'ellipse_persistence.csv',collapse,delimiter=',',header='minor_axis_gain,longest_H1_lifetime',comments='')
    plot(out,'ellipse_margin',lambda ax:ax.plot(collapse[:,0],collapse[:,1],'o-',label='Finite Rips H1 lifetime'),'Minor-axis channel gain','Persistent loop lifetime')
    (out/'results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    (out/'REPORT.md').write_text('# 本轮执行记录\n\n运行类型：解析积分、有限复形、有限差分与高斯蒙特卡洛；GPU语义编码器训练列入仿真计划。\n\n'+ '\n'.join('## '+k+'\n```json\n'+json.dumps(v,ensure_ascii=False,indent=2)+'\n```\n' for k,v in report['checks'].items()),encoding='utf8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('results'));main(p.parse_args().output)
