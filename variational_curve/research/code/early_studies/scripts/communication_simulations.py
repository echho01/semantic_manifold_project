#!/usr/bin/env python3
"""CPU numerical studies for task-oriented geometric/topological communication.
All curves are executed calculations, not image-model or GPU benchmark results.
Run from repo: OPENBLAS_NUM_THREADS=1 python scripts/communication_simulations.py
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1'); os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, csv, json, math, platform, time, itertools
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar, linear_sum_assignment
from scipy.special import logsumexp, xlogy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LN2=np.log(2.)

def savecsv(path, rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def figsave(root,name):
    ax=plt.gca();ax.grid(True,alpha=.22);ax.tick_params(labelsize=12)
    plt.tight_layout(pad=.7);plt.savefig(root/(name+'.pdf'),bbox_inches='tight');plt.savefig(root/(name+'.png'),dpi=220,bbox_inches='tight');plt.close()

def sphere(m,s):
    # Sphere radius one; theta has density proportional sin(theta)^(m-1).
    vmax=min(np.pi*np.sqrt(s),12.) if s>0 else np.pi
    if s>0:
        den=quad(lambda t: np.exp(-t*t)*np.sin(t/np.sqrt(s))**(m-1)/np.sqrt(s),0,vmax,epsabs=1e-14,epsrel=2e-11)[0]
        num=quad(lambda t: t*t/s*np.exp(-t*t)*np.sin(t/np.sqrt(s))**(m-1)/np.sqrt(s),0,vmax,epsabs=1e-14,epsrel=2e-11)[0]
    else:
        den=quad(lambda t: np.sin(t)**(m-1),0,np.pi)[0]
        num=quad(lambda t: t*t*np.sin(t)**(m-1),0,np.pi)[0]
    den0=quad(lambda t: np.sin(t)**(m-1),0,np.pi)[0]
    d=num/den;r=(-np.log(den/den0)-s*d)/LN2
    return d,max(0.,r)

def ba(w,c,lam,tol=2e-10,maxiter=10000):
    # Natural-log BA; finite source/reproduction alphabets, not a continuous RDF.
    pi=np.ones(c.shape[1])/c.shape[1]
    for it in range(maxiter):
        lq=np.log(np.maximum(pi,1e-300))[None,:]-lam*c
        lq-=logsumexp(lq,axis=1)[:,None];q=np.exp(lq)
        new=w@q
        if np.max(np.abs(new-pi))<tol: pi=new;break
        pi=new
    lq=np.log(np.maximum(pi,1e-300))[None,:]-lam*c
    lz=logsumexp(lq,axis=1);q=np.exp(lq-lz[:,None]);piout=w@q
    d=float(np.sum(w[:,None]*q*c));info=float(np.sum(w[:,None]*q*(np.log(np.maximum(q,1e-300))-np.log(np.maximum(piout,1e-300))[None,:])))
    # First-order dual gap for the mixture optimization: max_j sum_i w_i Kij/Zi - 1.
    dual=float(max(0.,np.max(np.exp(-lam*c-lz[:,None]).T@w)-1.))
    return info,d,it+1,dual

def geometry(root):
    rows=[]
    plt.figure(figsize=(6.3,4.2))
    for m in [1,2,4]:
        d0,_=sphere(m,0.)
        series=[]
        for s in np.geomspace(.04,150.,85):
            d,r=sphere(m,s);series.append((d/d0,r));rows.append(dict(m=m,s=float(s),D=d,D_normalized=d/d0,R_bits=r))
        a=np.array(series);plt.plot(a[:,0],a[:,1],lw=2,label=f'Intrinsic dimension m={m}')
    plt.xscale('log');plt.xlabel('Normalized semantic distortion D / D(0)',fontsize=13);plt.ylabel('Required rate (bit/state)',fontsize=13);plt.legend(fontsize=11);figsave(root,'G1_dimension_rate');savecsv(root/'G1_dimension_rate.csv',rows)
    # Exact sphere vs leading approximation; error of leading high-resolution formula.
    m=2;d0,_=sphere(2,0.);r2=[]
    for s in np.geomspace(.6,1000,90):
        d,r=sphere(2,s);b=(np.log(4*np.pi)-np.log(np.pi*np.e*d))/LN2
        r2.append(dict(s=float(s),D=d,R_exact=r,R_leading=b,R_curvature=b+d/(6*LN2)))
    plt.figure(figsize=(6.3,4.2));a=np.array([[r['D'],r['R_exact'],r['R_leading'],r['R_curvature']] for r in r2]);
    plt.plot(a[:,0],a[:,1],lw=2,label='Exact homogeneous-source RDF');plt.plot(a[:,0],a[:,2],'--',lw=1.8,label='Leading high-resolution term');plt.plot(a[:,0],a[:,3],':',lw=2,label='Including curvature correction');plt.xscale('log');plt.xlabel('Squared geodesic distortion D',fontsize=13);plt.ylabel('Required rate (bit/state)',fontsize=13);plt.legend(fontsize=10);figsave(root,'G2_curvature_rate');savecsv(root/'G2_curvature_rate.csv',r2)
    # A trainable, four-coefficient nonuniform circle parameterization, fixed semantic grid.
    n=96;x=2*np.pi*(np.arange(n)+.5)/n;w=(1+.85*np.cos(x))/n;dist=np.abs(x[:,None]-x[None,:]);cg=np.minimum(dist,2*np.pi-dist)**2
    basis=np.stack([np.cos(x),np.cos(2*x),np.sin(x),np.sin(2*x)],axis=1)
    def kernel(par):
        speed=np.exp(basis@par[:4]);edge=2*np.pi*speed/speed.sum();pos=np.cumsum(edge)-.5*edge
        dh=np.abs(pos[:,None]-pos[None,:]);ch=np.minimum(dh,2*np.pi-dh)**2
        s=np.exp(par[4]);lq=np.log(edge)[None,:]-s*ch;lq-=logsumexp(lq,axis=1)[:,None];q=np.exp(lq);pi=w@q
        d=np.sum(w[:,None]*q*cg);info=np.sum(w[:,None]*q*(lq-np.log(np.maximum(pi,1e-300))[None,:]))
        return float(info),float(d)
    out=[];history=[]
    for lam in np.geomspace(.25,32,25):
        scalar=minimize_scalar(lambda ls: sum(a*b for a,b in zip(kernel(np.r_[np.zeros(4),ls]),[1.,lam])),bounds=(-9,8),method='bounded',options={'xatol':1e-9})
        ini=np.r_[np.zeros(4),scalar.x]
        def obj(par):
            i,d=kernel(par);return i+lam*d
        trace=[obj(ini)]
        opt=minimize(obj,ini,method='L-BFGS-B',bounds=[(-.8,.8)]*4+[(-9,8)],callback=lambda p:trace.append(obj(p)),options={'maxiter':180,'ftol':1e-12,'gtol':2e-7})
        par=opt.x if opt.fun<=obj(ini) else ini
        fi,fd=kernel(ini);oi,od=kernel(par);bi,bd,nit,gap=ba(w,cg,lam)
        out.append(dict(lam=float(lam),fixed_R=fi/LN2,fixed_D=fd,optimized_R=oi/LN2,optimized_D=od,BA_R=bi/LN2,BA_D=bd,BA_dual_gap=gap,objective_fixed=obj(ini),objective_opt=obj(par),converged=bool(opt.success),a0=float(par[0]),a1=float(par[1]),a2=float(par[2]),a3=float(par[3]),log_s=float(par[4])))
        if len(history)==0 or abs(np.log(lam/2.))<abs(np.log(history[0]['lam']/2.)):
            history=[dict(lam=float(lam),iteration=i,objective=float(o)) for i,o in enumerate(trace)]
    plt.figure(figsize=(6.3,4.2))
    for prefix,label,style in [('fixed','Fixed geometry Gibbs code','--'),('optimized','Variational geometry Gibbs code','-'),('BA','Finite-alphabet BA benchmark',':')]:
        a=np.array([[r[prefix+'_D'],r[prefix+'_R']] for r in out]);plt.plot(a[:,0],a[:,1],style,lw=2,label=label)
    plt.xlabel('Fixed semantic distortion D',fontsize=13);plt.ylabel('Achievable rate (bit/state)',fontsize=13);plt.xlim(0,1.1);plt.legend(fontsize=10);figsave(root,'G3_variational_coding');savecsv(root/'G3_variational_coding.csv',out);savecsv(root/'G4_optimization.csv',history)
    plt.figure(figsize=(6.3,4.2));plt.plot([h['iteration'] for h in history],[h['objective'] for h in history],'-o',markersize=3);plt.xlabel('Optimization iteration',fontsize=13);plt.ylabel('I + lambda D (nat/state)',fontsize=13);figsave(root,'G4_optimization')
    # interpolate all curves only inside their shared domain; report matched-D gain.
    ds=.30
    def interpolate(prefix):
        pairs=sorted([(r[prefix+'_D'],r[prefix+'_R']) for r in out]);return float(np.interp(ds,*np.array(pairs).T))
    return {'dimension':{str(m):{'D0':sphere(m,0)[0],'R_at_normalized_D_0.1':float(np.interp(.1,*np.array(sorted([(r['D_normalized'],r['R_bits']) for r in rows if r['m']==m])).T))} for m in [1,2,4]},'geometry_coding':{'n':n,'source':'p(theta)=(1+0.85 cos theta)/(2pi)','parameter_bounds':[-.8,.8],'matched_D':ds,'fixed_bits':interpolate('fixed'),'optimized_bits':interpolate('optimized'),'BA_bits':interpolate('BA'),'max_BA_dual_gap':max(r['BA_dual_gap'] for r in out),'optimizer_successes':sum(r['converged'] for r in out),'operational_scope':'Explicit test-channel rates; not a measured entropy-coded bitstream.'}}

def hb(p):
    p=np.clip(p,0.,1.);return float(-(xlogy(p,p)+xlogy(1-p,1-p))/LN2)

def knowledge(root):
    rows=[];plt.figure(figsize=(6.3,4.2))
    for rho in [0.,.5,.8]:
        a=[]
        for eps in np.linspace(0,.5,101):
            val=(1-rho)*(1-hb(min(.5,eps/(1-rho))))
            rows.append(dict(known_probability=rho,structural_error_budget=float(eps),additional_bits=val));a.append(val)
        plt.plot(np.linspace(0,.5,101),a,lw=2,label=f'Shared type knowledge probability = {rho:g}')
    plt.xlabel('Allowed structural-type error probability',fontsize=13);plt.ylabel('Additional type rate (bit/episode)',fontsize=13);plt.legend(fontsize=10);figsave(root,'T1_structure_rate');savecsv(root/'T1_structure_rate.csv',rows)
    # independent binary type and circular coordinate: joint BA verifies tensorization.
    n=32;ang=2*np.pi*np.arange(n)/n;dd=np.abs(ang[:,None]-ang[None,:]);dc=np.minimum(dd,2*np.pi-dd)**2
    typ=np.repeat([0,1],n);fullc=np.tile(dc,(2,2));wc=np.ones(n)/n;wf=np.ones(2*n)/(2*n)
    checks=[]
    for l in [.5,2,6]:
        for b in [.4,1.5,4.]:
            i,d,_,gap=ba(wf,l*fullc+b*(typ[:,None]!=typ[None,:]),1.)
            ic,ddc,_,_=ba(wc,dc,l);err=1/(1+np.exp(b));expected=ic/LN2+1-hb(err)
            checks.append(dict(lambda_g=l,lambda_t=b,joint_BA_bits=i/LN2,product_bits=expected,abs_error=abs(i/LN2-expected),dual_gap=gap))
    savecsv(root/'T2_joint_BA_check.csv',checks)
    # Counterexample: global permutation preserves marginal and distances but swaps all labels.
    mis=[]
    for angle in np.linspace(0,np.pi,13):
        x=2*np.pi*(np.arange(4096)+.5)/4096;y=(x+angle)%(2*np.pi)
        pe=np.mean((x<np.pi)!=(y<np.pi));mis.append(dict(rotation_rad=float(angle),task_error=float(pe),marginal_TV=0.,diagram_distance=0.))
    savecsv(root/'T3_mismatch_counterexample.csv',mis)
    return {'type_rate_example':{'rho':.5,'epsilon':.05,'bits':.5*(1-hb(.1))},'tensorization_max_error':max(x['abs_error'] for x in checks),'mismatch_counterexample':'Uniform circle rotated by pi: identical distribution/PH, error one for a fixed binary semantic partition.'}

def rips(points):
    n=len(points);dist=np.linalg.norm(points[:,None]-points[None,:],axis=-1)
    cells=[(0.,(i,)) for i in range(n)]
    for size in (2,3):
        for s in itertools.combinations(range(n),size):cells.append((max(dist[i,j] for i,j in itertools.combinations(s,2)),s))
    cells.sort(key=lambda v:(v[0],len(v[1]),v[1]));index={s:i for i,(_,s) in enumerate(cells)};piv={};diags=[[],[]]
    for death,s in cells:
        col=0
        if len(s)>1:
            for f in itertools.combinations(s,len(s)-1):col ^= 1<<index[f]
        while col:
            k=col.bit_length()-1
            if k in piv:col^=piv[k]
            else:break
        if col:
            k=col.bit_length()-1;piv[k]=col;birth,bs=cells[k];dim=len(bs)-1
            if dim<2 and death>birth+1e-10:diags[dim].append((birth,death))
    return [np.array(d).reshape(-1,2) for d in diags]

def bn(a,b):
    n,m=len(a),len(b)
    if n+m==0:return 0.
    c=np.full((n+m,n+m),np.inf)
    if n and m:c[:n,:m]=np.max(np.abs(a[:,None]-b[None,:]),axis=-1)
    for i in range(n):c[i,m+i]=(a[i,1]-a[i,0])/2
    for j in range(m):c[n+j,j]=(b[j,1]-b[j,0])/2
    c[n:,m:]=0.;v=np.unique(c[np.isfinite(c)]);lo=0;hi=len(v)-1
    while lo<hi:
        mid=(lo+hi)//2;r,k=linear_sum_assignment((c>v[mid]).astype(float))
        if np.all(c[r,k]<=v[mid]):hi=mid
        else:lo=mid+1
    return float(v[lo])

def phdist(a,b):return max(bn(a[0],b[0]),bn(a[1],b[1]))
def wilson(k,n):
    z=1.95996398454;p=k/n;den=1+z*z/n;ctr=(p+z*z/(2*n))/den;rad=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;return ctr-rad,ctr+rad

def wireless(root,trials=350,seed=20260916):
    rng=np.random.default_rng(seed);n=12;theta=2*np.pi*(np.arange(n)+.13)/n
    # Same second moments and power convention; only relational shape differs.
    ring=np.stack([np.cos(theta),np.sin(theta)],axis=1)
    t=2*np.pi*(np.arange(n//2)+.13)/(n//2)
    cl=np.concatenate([np.stack([-.85+.22*np.cos(t),.22*np.sin(t)],1),np.stack([.85+.22*np.cos(t),.22*np.sin(t)],1)])
    clouds=[]
    for z in [ring,cl]:
        z=z-z.mean(0);z=z/np.sqrt(np.mean(z*z,axis=0));clouds.append(z)
    prototypes=[rips(z) for z in clouds];margin=phdist(*prototypes)
    # Task head: nearest prototype; certificate tau=prototype separation/2.
    tau=margin/2;alpha=.05;talpha=np.log(n*(n-1)/(2*alpha));h=np.array([1.,.22]);c=1/h**2
    allocs={'equal':np.array([.5,.5]),'certificate':np.sqrt(c)/np.sqrt(c).sum()}
    rows=[];sample=[]
    for snr in np.arange(0,43,3):
        P=2*10**(snr/10) # mean transmitted energy per real coordinate / noise variance.
        errors={k:0 for k in allocs};bounds={};viol={k:0 for k in allocs};dbs={k:[] for k in allocs}
        noises=rng.normal(size=(trials,n,2));types=rng.integers(0,2,size=trials)
        for name,share in allocs.items():
            p=P*share;sigma=c/p;eta=np.sqrt(2*sigma.sum())+2*np.sqrt(sigma.max()*talpha);bounds[name]=eta
            for i in range(trials):
                z=clouds[types[i]];zh=z+noises[i]*np.sqrt(sigma)[None,:];dh=rips(zh)
                d0=phdist(dh,prototypes[0]);d1=phdist(dh,prototypes[1]);pred=int(d1<d0);errors[name]+=int(pred!=types[i]);db=d0 if types[i]==0 else d1;dbs[name].append(db);viol[name]+=int(db>eta)
                if snr in [9,18,27] and i<24:
                    dmat=np.linalg.norm(z[:,None]-z[None,:],axis=-1);dhat=np.linalg.norm(zh[:,None]-zh[None,:],axis=-1);pair=np.max(np.abs(dmat-dhat));sample.append(dict(snr_db=int(snr),method=name,trial=i,PH_distance=float(db),pair_bound=float(pair),prob_radius=float(eta)))
            lo,hi=wilson(errors[name],trials)
            rows.append(dict(snr_db=int(snr),method=name,P_total=float(P),error_count=errors[name],trials=trials,error_rate=errors[name]/trials,CI_low=lo,CI_high=hi,PH_mean=float(np.mean(dbs[name])),radius=float(eta),certified=int(eta<tau),radius_violation=viol[name]))
    plt.figure(figsize=(6.3,4.2))
    for name,label in [('equal','Equal power'),('certificate','Certificate-based power')]:
        rr=[r for r in rows if r['method']==name];x=np.array([r['snr_db'] for r in rr]);y=np.array([r['error_rate'] for r in rr]);lo=np.array([r['CI_low'] for r in rr]);hi=np.array([r['CI_high'] for r in rr]);plt.plot(x,y,'-o',lw=2,ms=3.5,label=label);plt.fill_between(x,lo,hi,alpha=.13)
    plt.xlabel('Transmit SNR (dB)',fontsize=13);plt.ylabel('Structural decision error probability',fontsize=13);plt.legend(fontsize=11);figsave(root,'T4_wireless_error');savecsv(root/'T4_wireless_error.csv',rows);savecsv(root/'T5_stability_samples.csv',sample)
    plt.figure(figsize=(6.3,4.2))
    for name,label in [('equal','Equal-power radius'),('certificate','Allocated-power radius')]:
        rr=[r for r in rows if r['method']==name];plt.semilogy([r['snr_db'] for r in rr],[r['radius'] for r in rr],lw=2,label=label)
    plt.axhline(tau,ls='--',lw=1.8,label='Decision half-margin');plt.xlabel('Transmit SNR (dB)',fontsize=13);plt.ylabel('High-probability PH error radius',fontsize=13);plt.legend(fontsize=10);figsave(root,'T6_wireless_certificate')
    # Certified sufficient transmit power vs fixed worst-mode gain, two allocation schemes.
    powers=[];K=1+np.sqrt(2*talpha);slack=.01*tau;T0=(tau-slack)**2/(2*K*K)
    for hw in np.geomspace(.05,1,60):
        cc=np.array([1.,1/hw**2]);pop=np.sqrt(cc).sum()**2/T0;peq=2*cc.sum()/T0;powers.append(dict(weak_gain=float(hw),P_equal=float(peq),P_certificate=float(pop),gain_dB=float(10*np.log10(peq/pop))))
    plt.figure(figsize=(6.3,4.2));plt.plot([p['weak_gain'] for p in powers],[10*np.log10(p['P_equal']/2) for p in powers],'--',lw=2,label='Equal power');plt.plot([p['weak_gain'] for p in powers],[10*np.log10(p['P_certificate']/2) for p in powers],lw=2,label='Minimum for trace certificate');plt.xscale('log');plt.xlabel('Weak-mode amplitude gain',fontsize=13);plt.ylabel('Sufficient transmit SNR (dB)',fontsize=13);plt.legend(fontsize=10);figsave(root,'T7_required_power');savecsv(root/'T7_required_power.csv',powers)
    return {'trials_per_snr_method':trials,'N':n,'noise_variance':1.,'gains':h.tolist(),'power_shares':{k:v.tolist() for k,v in allocs.items()},'prototype_separation':margin,'decision_half_margin':tau,'alpha':alpha,'num_transmissions':len(rows)*trials,'max_pair_bound_violation':max(r['PH_distance']-r['pair_bound'] for r in sample),'total_radius_violations':sum(r['radius_violation'] for r in rows),'certification_SNR':{k:next((r['snr_db'] for r in rows if r['method']==k and r['certified']),None) for k in allocs},'selected_rows':[r for r in rows if r['snr_db'] in [9,12,15,18,21,24]],'scope':'Finite two-class synthetic episodes and known deterministic diagonal fading; no learned-image task or outage generalization.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/communication');p.add_argument('--trials',type=int,default=350);args=p.parse_args();root=Path(args.out);root.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.labelsize':13,'axes.titlesize':14,'legend.frameon':False})
    t=time.time();report={'kind':'Executed CPU analytic and synthetic-codec simulations','seed':20260916,'python':platform.python_version(),'numpy':np.__version__,'trials':args.trials}
    report.update(geometry(root));print('geometry completed',flush=True);report.update(knowledge(root));print('knowledge completed',flush=True);report['wireless']=wireless(root,args.trials);report['elapsed_seconds']=time.time()-t
    (root/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
