#!/usr/bin/env python3
"""Executed refinement: distinguish geometry deformation from codeword density.
The test channel changes; the source and semantic distortion NEVER change.
Rates are mutual informations, not measured bitstream lengths.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize
from communication_simulations import savecsv,figsave,ba,LN2
import matplotlib.pyplot as plt

def main(out):
    torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
    root=Path(out);root.mkdir(parents=True,exist_ok=True)
    old=np.genfromtxt(root/'G3_variational_coding.csv',names=True,delimiter=',',dtype=None,encoding=None)
    n=96;x=2*np.pi*(np.arange(n)+.5)/n
    w=torch.tensor((1+.85*np.cos(x))/n)
    dif=np.abs(x[:,None]-x[None,:]);cg=torch.tensor(np.minimum(dif,2*np.pi-dif)**2)
    bg=torch.tensor(np.stack([np.cos(x),np.cos(2*x),np.sin(x),np.sin(2*x)],1))
    br=torch.tensor(np.stack([fun(k*x) for k in range(1,5) for fun in (np.cos,np.sin)],1))
    def calc(par,geo):
        if geo:
            a=par[:4];b=par[4:12];ls=par[-1]
        else:a=par[:4]*0;b=par[4:12];ls=par[-1]
        speed=torch.exp(bg@a);edge=2*np.pi*speed/speed.sum();pos=torch.cumsum(edge,0)-.5*edge
        d=torch.abs(pos[:,None]-pos[None,:]);cost=torch.minimum(d,2*np.pi-d)**2
        lq=torch.log(edge)[None,:]+(br@b)[None,:]-torch.exp(ls)*cost
        lq=lq-torch.logsumexp(lq,1)[:,None];q=torch.exp(lq);pi=w@q
        I=torch.sum(w[:,None]*q*(lq-torch.log(pi.clamp_min(1e-300))[None,:]));D=torch.sum(w[:,None]*q*cg)
        return I,D
    rows=[];traceout=[];start=time.time()
    for row in old:
        lam=float(row['lam']);r={name:row[name].item() if hasattr(row[name],'item') else row[name] for name in old.dtype.names}
        for geo,key in [(False,'density'),(True,'joint')]:
            init=np.zeros(13);init[-1]=float(row['log_s']);init[:4]=[row[f'a{k}'] for k in range(4)] if geo else 0
            history=[]
            def obj(par):
                p=torch.tensor(par,requires_grad=True);i,d=calc(p,geo);v=i+lam*d;v.backward();return float(v.detach()),p.grad.numpy()
            def cb(p):history.append(obj(p)[0])
            op=minimize(obj,init,jac=True,method='L-BFGS-B',bounds=[(-.8,.8)]*4+[(-5,5)]*8+[(-9,8)],callback=cb,options={'maxiter':500,'ftol':2e-13,'gtol':1e-7})
            with torch.no_grad():i,d=calc(torch.tensor(op.x),geo)
            r[key+'_R']=float(i)/LN2;r[key+'_D']=float(d);r[key+'_objective']=float(i+lam*d);r[key+'_success']=bool(op.success)
            if geo and abs(lam-2)<.25:traceout=[dict(iteration=k,objective=v,lam=lam) for k,v in enumerate(history)]
        rows.append(r)
    savecsv(root/'G5_geometry_density.csv',rows)
    if traceout:savecsv(root/'G5_objective.csv',traceout)
    plt.figure(figsize=(6.3,4.2))
    for prefix,label,style in [('fixed','Fixed geometry + uniform measure','--'),('optimized','Geometry only',':'),('joint','Geometry + coding density','-'),('BA','Finite-source BA benchmark','-.')]:
        plt.plot([r[prefix+'_D'] for r in rows],[r[prefix+'_R'] for r in rows],style,lw=2,label=label)
    plt.xlabel('Fixed semantic distortion D',fontsize=13);plt.ylabel('Achievable rate (bit/state)',fontsize=13);plt.xlim(.06,.9);plt.ylim(.35,2.8);plt.legend(fontsize=9.5);figsave(root,'G5_geometry_density')
    D=.3;matched={}
    for k in ['fixed','optimized','density','joint','BA']:
        pts=np.array(sorted((r[k+'_D'],r[k+'_R']) for r in rows));matched[k]=float(np.interp(D,pts[:,0],pts[:,1]))
    summary={'kind':'Executed finite-source test-channel optimization','N':n,'source_amplitude':.85,'fixed_semantic_distortion':D,'matched_rate_bits_interpolated':matched,'geometry_parameters':4,'log_density_parameters':8,'dual_multiplier_count':len(rows),'elapsed_seconds':time.time()-start,'converged_joint':sum(r['joint_success'] for r in rows),'converged_density':sum(r['density_success'] for r in rows),'interpretation':'Geometry-only gain is small; density accounts for most reduction. Do not attribute joint gain to curvature or report as neural codec performance.'}
    (root/'refinement_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/communication');main(p.parse_args().out)
