"""Frozen-background dynamical Hartree gap elastic pilot.

Recompute mass from the evolving represented gas. Each macro step solves
with an implicit conservative fixed-gap kernel. Reject macro steps when
the gap or stationary Hartree functional changes too much. No energy
projection, no condensate, and no assertion of tree kinetic validity
after lambda*f exceeds unity.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from elastic_implicit_trial import ImplicitElastic,merge_grid
from quartic_consistency import shell_nodes,lambda_from_sigma,gap,variance,MPL
from hartree_condensate_limits import initial_U

def evolve_dynamic(e,f,m0,lam,TH,end,dm_limit=.025,rtol=2e-5,wall=300.,validity_limit=1.):
    start=time.monotonic();vol=np.diff(e**3)/(6*np.pi**2)
    f=np.maximum(f,1e-38);N=float(vol@f);U0=initial_U(e,f,m0,lam)['U']
    q2mom=np.diff(e**5)/(10*np.pi**2)
    init=f.copy();g=initial_U(e,f,m0,lam);t=0.;dt=min(end,2e-7);hist=[];rej=0
    def record(t,f,g,**kw):
        return dict(t_H=t,m_eff=g['m_eff'],energy_relative_error=g['U']/U0-1,
           number_relative_error=float(vol@f/N-1),lambda_fmax=float(lam*max(f)),
           fraction_lambda_f_gt1=float(vol[f*lam>1]@f[f*lam>1]/N),
           q2_relative_change=float(q2mom@f/(q2mom@init)-1),
           low_two_fraction=float(vol[:2]@f[:2]/N),elapsed=time.monotonic()-start,**kw)
    hist.append(record(t,f,g))
    message='finished';scale=lam**2*TH/(2048*np.pi**5)
    while t<end*(1-1e-12):
        if time.monotonic()-start>wall:message='wall budget';break
        if g['lambda_f_peak']>=validity_limit:message='tree kinetic occupation guard';break
        dt=min(dt,end-t);m=g['m_eff'];model=ImplicitElastic(e,m);norm=np.maximum(f,1.)
        def rhs(s,y):return model.rates_jac(y*norm,False)*scale/norm
        def jac(s,y):return model.rates_jac(y*norm,True)*scale*norm[None,:]/norm[:,None]
        sol=solve_ivp(rhs,(0,dt),f/norm,method='BDF',rtol=rtol,atol=1e-12,jac=jac,max_step=dt/4)
        if not sol.success:message='implicit step failure: '+sol.message;break
        ff=sol.y[:,-1]*norm
        neg=float(vol@np.maximum(-ff,0)/N)
        gg=initial_U(e,ff,m0,lam)
        dm=abs(gg['m_eff']/m-1)
        de=abs(gg['U']/g['U']-1)
        # Both local conservation and gap step checked; accumulated error is reported.
        if dm>dm_limit or de>2e-5 or neg>1e-9:
            dt*=.5;rej+=1
            if dt<1e-17:message='macro-step underflow';break
            continue
        t+=dt;f=ff;g=gg
        hist.append(record(t,f,g,dt_H=dt,relative_gap_step=dm,local_energy_drift=de,nfev=sol.nfev))
        print('state',json.dumps(hist[-1]),flush=True)
        if dm<dm_limit/3 and de<5e-6:dt*=1.8
    return dict(success=message=='finished',message=message,requested_t_H=end,achieved_t_H=t,
                dm_limit=dm_limit,rejected_macro_steps=rej,final=hist[-1],history=hist,edges=e.tolist(),
                initial_f=init.tolist(),final_f=f.tolist())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--sigma',type=float,default=1e-7);p.add_argument('--n',type=int,default=64)
    p.add_argument('--qmin',type=float,default=1e-5);p.add_argument('--end',type=float,default=1e-4)
    p.add_argument('--dm-limit',type=float,default=.025);p.add_argument('--validity-limit',type=float,default=1.)
    p.add_argument('--wall',type=float,default=280.);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];h=row['final'];x=h['x']
    e,f=merge_grid(row['frozen_cell_spectrum']['edges'],row['frozen_cell_spectrum']['f_phi'],a.n,a.qmin)
    lam=float(lambda_from_sigma(a.sigma));TH=MPL*x/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV'])*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
    r=evolve_dynamic(e,f,1e-12*x,lam,TH,a.end,a.dm_limit,wall=a.wall,validity_limit=a.validity_limit)
    r.update(sigma=a.sigma,n=len(f),qmin=a.qmin,validity_limit=a.validity_limit)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print('FINAL',json.dumps(r['final']),r['message'],flush=True)
