"""Post-process an unshifted coupled trajectory with shifted kinematic apertures."""
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre, expit, spence

from early_pole_frozen_dynamic_vlasov_reference import load
from early_pole_separate_fermion_heat_reference import Separate
from early_pole_decay_mode_clock_reference import intervals,Q


def calc(d,g,ex,order):
    s=Separate(d,g,1.,order,ex,True)
    span=np.log(.052/s.x0)
    tol=np.r_[np.full(s.ns,3e-8),np.full(5,6e-13)]
    sol=solve_ivp(s.local,(0.,span),s.y0,method='BDF',jac=s.jac,atol=tol,rtol=2e-9,
                  max_step=.003,dense_output=True)
    assert sol.success,sol.message
    nodes,w=roots_legendre(20)
    names=('bare','common_momentum_fixed','flavor_momentum_fixed','flavor_fixed_energy_control')

    def sample(t):
        x=s.x0*np.exp(t)
        a=sol.sol(t)
        m=s.d.d.gap(x,a[:s.ns])
        nb,np_,rb,rp,_=a[-5:]
        _,ab,tb=s.fermion(nb,rb,s.d.daughter_energies(x),s.d.d.wf)
        _,ap,tp=s.fermion(np_,rp,s.d.parent_energies(x),s.d.pw)
        shifts=1/6-np.array([tp*tp*spence(1+np.exp(ap)),tb*tb*spence(1+np.exp(ab))])/np.pi**2
        shifts*=g*g
        vals={}
        for label,dp,db in zip(names,(0.,shifts.mean(),shifts[0],shifts[0]),
                               (0.,shifts.mean(),shifts[1],shifts[1]),strict=True):
            es,lo,hi,_,_=intervals(np.hypot(x,np.sqrt(dp)),m,
                                  np.hypot(s.d.daughter_mass(x),np.sqrt(db)),Q)
            P=(lo+hi)/2+(hi-lo)*nodes/2
            B=P-es
            if label!='flavor_fixed_energy_control':
                P=np.sqrt(P*P-dp)
                B=np.sqrt(B*B-db)
            vals[label]=float(.5*w@(expit(ap-P/tp)-expit(ab-B/tb)))
        return dict(x=x,photon_moment_fixed=1/6,mean_scalar_mass=m,temperature=[tb,tp],
                    alpha=[ab,ap],difference_mass_sq=shifts[0]-shifts[1],gate=vals)

    t=np.linspace(0.,span,111)
    rows=[sample(v) for v in t]
    roots={}
    for name in names:
        iii=np.flatnonzero([r['gate'][name]<=0. for r in rows])
        i=iii[0] if len(iii) else None
        roots[name]=None if i is None else float(s.x0*np.exp(
            brentq(lambda v:sample(v)['gate'][name],t[i-1],t[i],xtol=1e-13))) if i else s.x0
    return dict(g=g,exchange=ex,order=order,scope='Post-processing only; the entire '
                'kinetic feedback and fitted moments are the unshifted Separate model.',
                roots=roots,rows=rows,final_scalar_count=float(s.weights@sol.y[:s.ns,-1]))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--g',type=float,default=.005)
    p.add_argument('--exchange',type=float,default=100.)
    p.add_argument('--order',type=int,default=10)
    o=p.parse_args()
    d=load('results/early-pole-spectrum-third-m06-2026-11-26.json.gz')
    print(json.dumps(calc(d,o.g,o.exchange,o.order),indent=2,allow_nan=False))
