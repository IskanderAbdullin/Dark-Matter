"""Necessary parameter tradeoff for closing inelastic halo channels.

Use the soft-mixing EFT of the preceding checkpoint. Impose illustrative
radiative splitting floor, log vacuum threshold and thermal/vacuum energy
cap. No claim of a regulator-independent naturalness theorem. Optimize
over full light-state B fraction and permitted Delta, not just UOPT.
Maxwell pair rates are exact for the stated NR equal-mass contact model.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.special import kve,zeta
from scipy.integrate import quad
from quartic_consistency import lambda_from_sigma
from direct_activation_consistency import radiation_T4,rho_tot
from resonant_mediator_inverse import M,C

def threshold_factor(d):
    if d<1e-10:return 1.
    return .5*d*np.exp(-d)*kve(1,.5*d)
def ratios(u,d):
    r=(1-u)/u
    return 2*r*threshold_factor(d),r*r*threshold_factor(2*d)
def sumratio(u,d):return sum(ratios(u,d))
def gcrit(d):
    if d<=1e-12:return 2.
    return 1+quad(lambda t:2*t*t*np.exp(-t*t-d)/(-np.expm1(-t*t-d)),0,10,epsabs=1e-13,epsrel=1e-9)[0]/(np.sqrt(np.pi)/2*zeta(1.5,1))
def min_delta(u,lam):return lam*M*M/(32*np.pi**2*u*u)
def loop_X(u,d,T):
    delta=2*M*d*T*np.sqrt(u*(1-u))
    return delta**2/(16*np.pi**2)
def Xcap(z,dneff):
    return min(np.pi**2/9*radiation_T4(z,dneff),.04*rho_tot(z))
def optimize(sigma,V,z,dneff=.1):
    lam=float(lambda_from_sigma(sigma));T=M*(V*1e5/C)**2
    cap=Xcap(z,dneff);delta_max=4*np.pi*np.sqrt(cap)
    uopt=(3-np.sqrt(3))/4
    delta_floor=lambda u:min_delta(u,lam)*np.sqrt(u*(1-u))/(1-2*u)
    margin=delta_max/delta_floor(uopt)
    result=dict(sigma=sigma,V0_kms=V,z_on=z,dneff=dneff,T_eV=T,Xcap_ev4=cap,delta_max=delta_max,overlap_exists=bool(margin>=1),max_margin_delta=margin)
    if margin>=1:
        ulo=brentq(lambda u:delta_floor(u)/delta_max-1,1e-9,uopt)
        uhi=brentq(lambda u:delta_floor(u)/delta_max-1,uopt,.5-1e-9)
        dmax=lambda u:delta_max/(2*M*T*np.sqrt(u*(1-u)))
        U=np.linspace(ulo,uhi,1200)
        rate=np.array([sumratio(u,dmax(u)) for u in U])
        ix=int(np.argmin(rate))
        result.update(u_range=[ulo,uhi],maximum_gap_over_T=float(dmax(ulo)),
           minimum_total_conversion_over_LL=float(rate[ix]),min_rate_u=float(U[ix]),
           minimum_critical_density_ratio=float(gcrit(dmax(ulo))))
    requirements=[]
    d_crit=brentq(lambda d:gcrit(d)-1-.05,0.,100.)
    for target in [.1,.01]:
        def need(u):
            d_rate=brentq(lambda d:sumratio(u,d)-target,0.,100.)
            d_floor=min_delta(u,lam)/(2*M*T*(1-2*u))
            d=max(d_rate,d_crit,d_floor)
            return loop_X(u,d,T),d,d_rate,d_crit,d_floor
        # The min is nonsmooth at the active constraint; grid first,
        # then bounded refinement around the best grid point.
        U=np.linspace(.00001,.49999,180)
        a=np.array([need(u)[0] for u in U]);j=int(a.argmin())
        lo=U[max(j-1,0)];hi=U[min(j+1,len(U)-1)]
        opt=minimize_scalar(lambda u:need(u)[0],bounds=(lo,hi),method='bounded',options={'xatol':1e-11})
        X,d,dr,dc,df=need(opt.x)
        requirements.append(dict(target_conversion_over_LL=target,target_critical_density_shift=.05,
           optimized_u=float(opt.x),required_gap_over_T=float(d),required_gap_eV=float(d*T),
           required_X_ev4=float(X),required_X_over_cap=float(X/cap),
           gap_requirements=[float(dr),float(dc),float(df)]))
    result['requirements']=requirements
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    out=[]
    for s in [.1,.3,1.25]:
        for V in [10.,20.,44.1,100.]:
            for z in [3.,5.,8.,10.,20.]:
                for dneff in [.1,.3]:
                    out.append(optimize(s,V,z,dneff))
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in out:
        if r['V0_kms']==44.1 and r['dneff']==.1:
            print(r['sigma'],r['z_on'],r.get('maximum_gap_over_T'),r.get('minimum_total_conversion_over_LL'),r['requirements'][0]['required_X_over_cap'])
