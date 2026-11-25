"""Exact invariant gain coefficient for A -> B + massless eta.
Frozen isotropic parents, B=0; radiation is NOT assumed empty.
For constant snapshot, integrated gain across redshifting line is
2*pi^2 Gamma_A (n_A-n_B)/(H pstar^3) including time dilation.
Only linear, frozen-reservoir diagnostic; nonlinear backreaction required.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import quad,trapezoid
from scipy.special import logsumexp
from late_decay_bose_screen import mesh,vgrid
from late_activation_overlap import H
from resonant_mediator_inverse import M,HBAR,SEC_GYR,C,HBARC
from direct_activation_consistency import T0,RHO0_CGS,CGS_RHO_TO_EV4
from quartic_consistency import lambda_from_sigma

def eta_gain(k,r,ps,m,gamma,edges,f,p0,x=None):
    # Invariant C_eta with exact piecewise-constant f_A integral.
    x=np.log(k/ps) if x is None else np.asarray(x)
    r=np.asarray(r);k=np.asarray(k)
    qmin=m*abs(np.sinh(x))/(p0*r)
    lo=np.maximum(edges[:-1],qmin[...,None])*p0*r[...,None]
    hi=edges[1:]*p0*r[...,None]
    de=np.maximum(hi-lo,0)*(hi+lo)/(np.hypot(hi,m)+np.hypot(lo,m))
    return m*m*gamma/(2*ps*k*k)*np.sum(de*f,axis=-1)

def piecequad(fn,edges,order=12):
    from scipy.special import roots_legendre
    x,w=roots_legendre(order);a,b=edges[:-1],edges[1:]
    y=(a[:,None]+b[:,None])/2+(b-a)[:,None]*x/2
    return float(np.sum(fn(y)*w*(b-a)[:,None]/2))

def fixed_tau(r,ps,gamma,edges,f,p0,m=M,order=12):
    p=edges*p0*r
    t=np.arcsinh(p/m);N=f@np.diff(p**3)/(6*np.pi**2)
    theo=2*np.pi**2*gamma*N/(H(r-1)*ps**3)
    a=np.unique(np.r_[-t[::-1],t])
    val=piecequad(lambda x:eta_gain(ps*np.exp(x),r,ps,m,gamma,edges,f,p0,x)/H(r-1),a,order)
    return theo,val,N

def crossing_tau(rc,ps,gamma,edges,f,p0,m=M,r_start=1e8,r_end=1.,order=12):
    # Exact along a redshifting mode even at semi-relativistic parent p.
    y=2*edges*p0*rc/m
    lo=max(np.log(r_end/rc),-.5*np.log1p(y[-1]))
    hi=np.log(r_start/rc)
    if y[-1]<1:hi=min(hi,-.5*np.log1p(-y[-1]))
    if hi<=lo:return 0.
    cuts=np.unique(np.r_[lo,hi,0.,-.5*np.log1p(y),-.5*np.log1p(-y[y<1])])
    cuts=cuts[(cuts>=lo)&(cuts<=hi)]
    return piecequad(lambda x:eta_gain(ps*np.exp(x),rc*np.exp(x),ps,m,gamma,edges,f,p0,x)/H(rc*np.exp(x)-1),cuts,order)
def rates(rc,ps,lifetime,N,sigma):
    g=1/(lifetime*SEC_GYR);lam=float(lambda_from_sigma(sigma))
    density=RHO0_CGS*CGS_RHO_TO_EV4/M*rc**3
    Hn=HBAR*H(rc-1);g2=8*np.pi*M*M*g*HBAR/ps
    # Fraction-one B medium scales; NOT self-consistent mass dispersion.
    mu=lam*density/(8*M*M)
    loss=2*np.sqrt(3)/(12288*np.pi)*lam**4*density**3/M**8/Hn
    return dict(n_all_eV3=density,eta_vacuum_cubic_eV=np.sqrt(g2),eta_loop_mass_scale_eV=np.sqrt(g2)/(4*np.pi),
                meanfield_over_recoil_fullB=mu/ps,quartic_coherent_loss_per_H_fullB=loss,
                gas_parameter=density*(lam/(32*np.pi*M))**3,meanfield_over_mass=mu/M,
                fullB_contact_sound_speed_over_c=np.sqrt(mu/M))

def run(row,lifetime,sigma):
    p0=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']*row['final']['peak_bin_q']
    e,_,f,meta=mesh(row,160,20.,0.)
    gamma=1/(lifetime*SEC_GYR)
    out=[]
    for ps in [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,.1]:
        for r in [1e8,3e7,1e7,3e6,1e6,3e5,1e5,3e4,1e4,1e3,100.,10.,2.]:
            t,v,N=fixed_tau(r,ps,gamma,e,f,p0)
            assert abs(v/t-1)<5e-8,(v,t,r,ps)
            actual=crossing_tau(r,ps,gamma,e,f,p0)
            check=crossing_tau(r,ps,gamma,e,f,p0,order=24)
            assert abs(check/actual-1)<5e-7,(r,ps,check,actual)
            # Frozen f_eta=e^actual-1. Fraction of available parent number
            # in one log bin of emerging eta (local estimate, not evolution).
            logf=actual+np.log(-np.expm1(-actual)) if actual>1e-6 else np.log(np.expm1(actual))
            logfrac=logf+np.log(ps**3/(2*np.pi**2*N))
            out.append(dict(recoil_eV=ps,r=r,tau_exact=t,tau_numeric=v,tau_expanding_frozen=actual,
                relative_crossing_difference=actual/t-1,crossing_quadrature_error=check/actual-1,log10_f_eta_frozen=logf/np.log(10),
                log10_eta_per_log_k_over_nA_frozen=logfrac/np.log(10),**rates(r,ps,lifetime,N,sigma)))
    return dict(lifetime_Gyr=lifetime,sigma_cm2_g=sigma,p0=p0,mesh=meta,cases=out,
        scope='Exact invariant frozen line gain, initial B=0. Exponentially large frozen gain diagnoses failure of the dilute eta approximation, not actual photon production or a constraint on a coupled nonlinear solution.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--lifetime',type=float,default=10.);p.add_argument('--sigma',type=float,default=.3);a=p.parse_args()
    out=run(json.loads(Path(a.input).read_text())[0],a.lifetime,a.sigma)
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in out['cases']:
        if r['r'] in [1e5,1e3,10.]:print('recoil',r['recoil_eV'],'r',r['r'],'tau',r['tau_exact'],'mf/gap',r['meanfield_over_recoil_fullB'],'loss/H',r['quartic_coherent_loss_per_H_fullB'])
