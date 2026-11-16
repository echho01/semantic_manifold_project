#!/usr/bin/env python3
"""Reproducible CPU studies. Test-channel rates are NOT measured bitstream rates.
Run: OPENBLAS_NUM_THREADS=1 python scripts/finalize_numerical_study.py
No GPU experiment is implied. All costs use a fixed source/semantic distortion.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,csv,json,itertools,time,platform
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize,minimize_scalar,linear_sum_assignment
from scipy.special import logsumexp,xlogy
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
LN2=np.log(2.)

def savecsv(path,rows):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def figsave(root,name):
    plt.gca().grid(True,alpha=.22);plt.tight_layout(pad=.8)
    for ext in ['pdf','svg','png']:plt.savefig(root/(name+'.'+ext),dpi=240,bbox_inches='tight')
    plt.close()
def ba(w,c,lam,tol=2e-10,maxiter=10000):
    pi=np.ones(c.shape[1])/c.shape[1]
    for it in range(maxiter):
        lq=np.log(np.maximum(pi,1e-300))[None,:]-lam*c;lq-=logsumexp(lq,axis=1)[:,None];q=np.exp(lq);new=w@q
        if np.max(abs(new-pi))<tol:pi=new;break
        pi=new
    lq=np.log(np.maximum(pi,1e-300))[None,:]-lam*c;lz=logsumexp(lq,axis=1);q=np.exp(lq-lz[:,None]);pout=w@q
    d=float(np.sum(w[:,None]*q*c));info=float(np.sum(w[:,None]*q*(np.log(np.maximum(q,1e-300))-np.log(np.maximum(pout,1e-300))[None,:])))
    gap=max(0.,float(np.max(np.exp(-lam*c-lz[:,None]).T@w)-1.))
    return info,d,it+1,gap

def sphere(m,s):
    den0=quad(lambda t:np.sin(t)**(m-1),0,np.pi)[0]
    if s>0:
        vmax=min(np.pi*np.sqrt(s),12.)
        fun=lambda t:np.exp(-t*t)*np.sin(t/np.sqrt(s))**(m-1)/np.sqrt(s)
        den=quad(fun,0,vmax,epsabs=1e-14,epsrel=2e-11)[0];num=quad(lambda t:t*t/s*fun(t),0,vmax,epsabs=1e-14,epsrel=2e-11)[0]
    else:den=den0;num=quad(lambda t:t*t*np.sin(t)**(m-1),0,np.pi)[0]
    d=num/den;return d,max(0.,(-np.log(den/den0)-s*d)/LN2)

def geometry(root):
    rows=[];plt.figure(figsize=(7,4.7))
    for m in [1,2,4]:
        d0,_=sphere(m,0.);a=[]
        for s in np.geomspace(.04,150.,85):
            d,r=sphere(m,s);a.append((d/d0,r));rows.append(dict(m=m,s=float(s),D=d,D_normalized=d/d0,R_bits=r))
        a=np.array(a);plt.plot(a[:,0],a[:,1],lw=2.5,label=f'm = {m}')
    plt.xscale('log');plt.xlabel('Normalized semantic distortion');plt.ylabel('Required rate (bit/state)');plt.legend(title='Intrinsic dimension');figsave(root,'G1_dimension_rate');savecsv(root/'G1_dimension_rate.csv',rows)
    r2=[]
    for s in np.geomspace(.6,1000,90):
        d,r=sphere(2,s);b=(np.log(4*np.pi)-np.log(np.pi*np.e*d))/LN2;r2.append(dict(s=float(s),D=d,R_exact=r,R_leading=b,R_curvature=b+d/(6*LN2)))
    plt.figure(figsize=(7,4.7))
    for k,lab,sty in [('R_exact','Exact sphere RDF','-'),('R_leading','Leading term','--'),('R_curvature','With curvature',':')]:plt.plot([r['D'] for r in r2],[r[k] for r in r2],sty,lw=2.5,label=lab)
    plt.xscale('log');plt.xlabel('Squared geodesic distortion D');plt.ylabel('Required rate (bit/state)');plt.legend();figsave(root,'G2_curvature_rate');savecsv(root/'G2_curvature_rate.csv',r2)
    torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
    n=96;x=2*np.pi*(np.arange(n)+.5)/n;wn=(1+.85*np.cos(x))/n;w=torch.tensor(wn);dif=np.abs(x[:,None]-x[None,:]);cn=np.minimum(dif,2*np.pi-dif)**2;cg=torch.tensor(cn)
    bg=torch.tensor(np.stack([np.cos(x),np.cos(2*x),np.sin(x),np.sin(2*x)],1));br=torch.tensor(np.stack([fun(k*x) for k in range(1,5) for fun in (np.cos,np.sin)],1))
    def calc(par,mode):
        a=par[:4] if mode in ['geometry','joint'] else par[:4]*0;b=par[4:12] if mode in ['density','joint'] else par[4:12]*0
        speed=torch.exp(bg@a);edge=2*np.pi*speed/speed.sum();pos=torch.cumsum(edge,0)-.5*edge;ds=torch.abs(pos[:,None]-pos[None,:]);cost=torch.minimum(ds,2*np.pi-ds)**2
        lq=torch.log(edge)[None,:]+(br@b)[None,:]-torch.exp(par[-1])*cost;lq=lq-torch.logsumexp(lq,1)[:,None];q=torch.exp(lq);pi=w@q
        I=torch.sum(w[:,None]*q*(lq-torch.log(pi.clamp_min(1e-300))[None,:]));D=torch.sum(w[:,None]*q*cg)
        return I,D
    result=[];hist=[]
    for lam in np.geomspace(.25,32,25):
        def scalar(ls):
            z=torch.zeros(13);z[-1]=ls;i,d=calc(z,'fixed');return float(i+lam*d)
        sc=minimize_scalar(scalar,bounds=(-9,8),method='bounded',options={'xatol':1e-9});p0=np.zeros(13);p0[-1]=sc.x
        r={'lambda':float(lam)};store={}
        for mode in ['fixed','geometry','density','joint']:
            ini=p0.copy()
            if mode=='joint':ini[:4]=store['geometry'][:4];ini[-1]=store['geometry'][-1]
            def obj(p):
                par=torch.tensor(p,requires_grad=True);i,d=calc(par,mode);v=i+lam*d;v.backward();return float(v.detach()),par.grad.numpy()
            tr=[]
            if mode=='fixed':out=ini;success=True
            else:
                op=minimize(obj,ini,jac=True,method='L-BFGS-B',bounds=[(-.8,.8)]*4+[(-5,5)]*8+[(-9,8)],callback=lambda p:tr.append(obj(p)[0]),options={'maxiter':500,'ftol':2e-13,'gtol':1e-7});out=op.x;success=bool(op.success)
            store[mode]=out
            with torch.no_grad():i,d=calc(torch.tensor(out),mode)
            r[mode+'_R']=float(i)/LN2;r[mode+'_D']=float(d);r[mode+'_objective']=float(i+lam*d);r[mode+'_success']=success
            if mode=='joint' and abs(lam-2)<.25:hist=[dict(iteration=j,objective=v,lam=float(lam)) for j,v in enumerate(tr)]
        i,d,it,gap=ba(wn,cn,lam);r.update(BA_R=i/LN2,BA_D=d,BA_iterations=it,BA_dual_gap=gap);result.append(r)
    savecsv(root/'G5_geometry_density.csv',result);savecsv(root/'G5_objective.csv',hist)
    plt.figure(figsize=(7.6,4.8))
    for pre,lab,sty in [('fixed','Fixed kernel','--'),('geometry','Geometry only',':'),('density','Density only','-'),('joint','Joint optimization','-.'),('BA','Finite-source BA','--')]:
        plt.plot([r[pre+'_D'] for r in result],[r[pre+'_R'] for r in result],sty,lw=2.2,label=lab)
    plt.xlim(.06,.9);plt.ylim(.35,2.8);plt.xlabel('Fixed semantic distortion D');plt.ylabel('Achievable rate (bit/state)');plt.legend(fontsize=14);figsave(root,'G5_geometry_density')
    plt.figure(figsize=(7,4.7));plt.plot([r['iteration'] for r in hist],[r['objective'] for r in hist],'-o',ms=3,lw=2);plt.xlabel('Optimization iteration');plt.ylabel('I + lambda D (nat/state)');figsave(root,'G6_objective')
    matched={}
    for pre in ['fixed','geometry','density','joint','BA']:
        a=np.array(sorted((r[pre+'_D'],r[pre+'_R']) for r in result));matched[pre]=float(np.interp(.3,a[:,0],a[:,1]))
    dm={}
    for m in [1,2,4]:
        a=np.array(sorted((r['D_normalized'],r['R_bits']) for r in rows if r['m']==m));dm[str(m)]=float(np.interp(.1,a[:,0],a[:,1]))
    return dict(dimension_rates_at_normalized_D_01=dm,matched_D=.3,matched_test_channel_bits=matched,max_BA_dual_gap=max(r['BA_dual_gap'] for r in result),joint_optimizer_successes=sum(r['joint_success'] for r in result))

def hb(p):
    p=np.clip(p,0.,1.);return float(-(xlogy(p,p)+xlogy(1-p,1-p))/LN2)
def knowledge(root):
    rows=[];plt.figure(figsize=(7,4.7))
    for rho in [0.,.5,.8]:
        a=[]
        for eps in np.linspace(0,.5,101):
            val=(1-rho)*(1-hb(min(.5,eps/(1-rho))));rows.append(dict(known_probability=rho,error_budget=float(eps),additional_bits=val));a.append(val)
        plt.plot(np.linspace(0,.5,101),a,lw=2.5,label=f'rho = {rho:g}')
    plt.xlabel('Allowed structural-type error');plt.ylabel('Additional rate (bit/state)');plt.legend(title='Known-type probability');figsave(root,'T1_structure_rate');savecsv(root/'T1_structure_rate.csv',rows)
    n=32;ang=2*np.pi*np.arange(n)/n;dd=np.abs(ang[:,None]-ang[None,:]);dc=np.minimum(dd,2*np.pi-dd)**2;typ=np.repeat([0,1],n);fullc=np.tile(dc,(2,2));checks=[]
    for l in [.5,2,6]:
        for b in [.4,1.5,4.]:
            i,d,_,gap=ba(np.ones(2*n)/(2*n),l*fullc+b*(typ[:,None]!=typ[None,:]),1.);ic,_,_,_=ba(np.ones(n)/n,dc,l);expected=ic/LN2+1-hb(1/(1+np.exp(b)));checks.append(dict(lambda_g=l,lambda_t=b,joint_BA_bits=i/LN2,product_bits=expected,abs_error=abs(i/LN2-expected),dual_gap=gap))
    savecsv(root/'T2_product_check.csv',checks)
    mismatch=[]
    for ang in np.linspace(0,np.pi,13):
        x=2*np.pi*(np.arange(4096)+.5)/4096;y=(x+ang)%(2*np.pi);mismatch.append(dict(rotation_rad=float(ang),task_error=float(np.mean((x<np.pi)!=(y<np.pi))),marginal_TV=0.,PH_distance=0.))
    savecsv(root/'T3_mismatch.csv',mismatch)
    return dict(type_rate_rho05_eps005=.5*(1-hb(.1)),product_max_error=max(x['abs_error'] for x in checks))

def rips(z):
    n=len(z);dist=np.linalg.norm(z[:,None]-z[None,:],axis=-1);cells=[(0.,(i,)) for i in range(n)]
    for size in [2,3]:
        for s in itertools.combinations(range(n),size):cells.append((max(dist[i,j] for i,j in itertools.combinations(s,2)),s))
    cells.sort(key=lambda v:(v[0],len(v[1]),v[1]));index={s:i for i,(_,s) in enumerate(cells)};piv={};diags=[[],[]]
    for death,s in cells:
        col=0
        if len(s)>1:
            for f in itertools.combinations(s,len(s)-1):col ^= 1<<index[f]
        while col:
            k=col.bit_length()-1
            if k in piv:col ^= piv[k]
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
    rng=np.random.default_rng(seed);n=12;theta=2*np.pi*(np.arange(n)+.13)/n;ring=np.stack([np.cos(theta),np.sin(theta)],1);t=2*np.pi*(np.arange(n//2)+.13)/(n//2);cl=np.concatenate([np.stack([-.85+.22*np.cos(t),.22*np.sin(t)],1),np.stack([.85+.22*np.cos(t),.22*np.sin(t)],1)])
    clouds=[]
    for z in [ring,cl]:z=z-z.mean(0);z=z/np.sqrt(np.mean(z*z,axis=0));clouds.append(z)
    prototypes=[rips(z) for z in clouds];margin=phdist(*prototypes);tau=margin/2;alpha=.05;talpha=np.log(n*(n-1)/(2*alpha));h=np.array([1.,.22]);c=1/h**2;allocs={'equal':np.array([.5,.5]),'certificate':np.sqrt(c)/np.sqrt(c).sum()};rows=[];samples=[];paired=[]
    for snr in np.arange(0,43,3):
        P=2*10**(snr/10);noise=rng.normal(size=(trials,n,2));types=rng.integers(0,2,size=trials);errs={}
        for name,share in allocs.items():
            sigma=c/(P*share);eta=np.sqrt(2*sigma.sum())+2*np.sqrt(sigma.max()*talpha);ee=[];ds=[];viol=0
            for i in range(trials):
                z=clouds[types[i]];zh=z+noise[i]*np.sqrt(sigma)[None,:];dh=rips(zh);d0=phdist(dh,prototypes[0]);d1=phdist(dh,prototypes[1]);ee.append(int(int(d1<d0)!=types[i]));db=d0 if types[i]==0 else d1;ds.append(db);viol+=int(db>eta)
                if snr in [9,18,27] and i<24:
                    dm=np.linalg.norm(z[:,None]-z[None,:],axis=-1);dhm=np.linalg.norm(zh[:,None]-zh[None,:],axis=-1);samples.append(dict(snr_db=int(snr),method=name,trial=i,PH_distance=float(db),pair_bound=float(np.max(abs(dm-dhm))),prob_radius=float(eta)))
            errs[name]=np.array(ee);lo,hi=wilson(sum(ee),trials);rows.append(dict(snr_db=int(snr),method=name,P_total=float(P),error_count=int(sum(ee)),trials=trials,error_rate=float(np.mean(ee)),CI_low=lo,CI_high=hi,PH_mean=float(np.mean(ds)),radius=float(eta),certified=int(eta<tau),radius_violation=viol))
        for i in range(trials):paired.append(dict(snr_db=int(snr),trial=i,type=int(types[i]),equal_error=int(errs['equal'][i]),allocated_error=int(errs['certificate'][i])))
    savecsv(root/'T4_wireless_error.csv',rows);savecsv(root/'T4_paired_trials.csv',paired);savecsv(root/'T5_stability_samples.csv',samples)
    plt.figure(figsize=(7,4.7))
    for name,lab in [('equal','Equal power'),('certificate','Structure-protection allocation')]:
        rr=[r for r in rows if r['method']==name];x=[r['snr_db'] for r in rr];y=[r['error_rate'] for r in rr];plt.plot(x,y,'-o',lw=2.5,ms=4,label=lab);plt.fill_between(x,[r['CI_low'] for r in rr],[r['CI_high'] for r in rr],alpha=.13)
    plt.xlabel('Transmit SNR (dB)');plt.ylabel('Structural decision error');plt.legend(fontsize=13.5);figsave(root,'T4_wireless_error')
    plt.figure(figsize=(7,4.7))
    for name,lab in [('equal','Equal-power bound'),('certificate','Allocated-power bound')]:
        rr=[r for r in rows if r['method']==name];plt.semilogy([r['snr_db'] for r in rr],[r['radius'] for r in rr],lw=2.5,label=lab)
    plt.axhline(tau,ls='--',lw=2,label='Decision half-margin');plt.xlabel('Transmit SNR (dB)');plt.ylabel('High-probability error radius');plt.legend(fontsize=13);figsave(root,'T6_wireless_certificate')
    powers=[];K=1+np.sqrt(2*talpha);T0=(.99*tau)**2/(2*K*K)
    for hw in np.geomspace(.05,1,60):
        cc=np.array([1.,1/hw**2]);po=np.sqrt(cc).sum()**2/T0;pe=2*cc.sum()/T0;powers.append(dict(weak_gain=float(hw),P_equal=float(pe),P_allocated=float(po),saving_dB=float(10*np.log10(pe/po))))
    savecsv(root/'T7_required_power.csv',powers);plt.figure(figsize=(7,4.7))
    for key,lab,sty in [('P_equal','Equal power','--'),('P_allocated','Minimum for sufficient bound','-')]:plt.plot([p['weak_gain'] for p in powers],[10*np.log10(p[key]/2) for p in powers],sty,lw=2.5,label=lab)
    plt.xscale('log');plt.xlabel('Weak-mode amplitude gain');plt.ylabel('Sufficient transmit SNR (dB)');plt.legend(fontsize=13);figsave(root,'T7_required_power')
    return dict(trials_per_snr_method=trials,N=n,noise_variance=1.,gains=h.tolist(),power_shares={k:v.tolist() for k,v in allocs.items()},prototype_separation=margin,decision_half_margin=tau,alpha=alpha,transmissions=len(rows)*trials,max_pair_bound_violation=max(r['PH_distance']-r['pair_bound'] for r in samples),radius_violations=sum(r['radius_violation'] for r in rows),certification_SNR={k:next((r['snr_db'] for r in rows if r['method']==k and r['certified']),None) for k in allocs},selected_rows=[r for r in rows if r['snr_db']==24])

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/communication_final');p.add_argument('--trials',type=int,default=350);args=p.parse_args();root=Path(args.out);root.mkdir(parents=True,exist_ok=True)
    if args.trials<=0:raise ValueError('--trials must be positive')
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.labelsize':17,'xtick.labelsize':14,'ytick.labelsize':14,'legend.fontsize':15,'legend.title_fontsize':15,'legend.frameon':False})
    start=time.time();report={'kind':'Executed CPU analytic and synthetic study; not GPU/image benchmarks','seed':20260916,'python':platform.python_version(),'numpy':np.__version__,'torch':torch.__version__}
    report['geometry']=geometry(root);print('Geometry complete',flush=True);report['knowledge']=knowledge(root);report['wireless']=wireless(root,args.trials);report['elapsed_seconds']=time.time()-start
    (root/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
