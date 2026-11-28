"""Test nonadditivity of a mean-pressure surrogate in its own controlled model."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp


def run(eos,files):
    b=Background(eos)
    Hm=np.sqrt(b.om)/b.D0
    rows=[np.loadtxt(p).T for p in files]
    num=np.asarray([np.trapezoid(F*q*q,q) for q,F in rows])
    q2=np.asarray([np.trapezoid(F*q**4,q) for q,F in rows])/num
    assert abs(np.average(q2[1:],weights=num[1:])/q2[0]-1.)<3e-7
    beta=np.r_[q2*(5/9)*(b.T0/mp)**2,0.]
    f=num[2]/num[1:].sum()
    k=np.array([.3,1.,2.,4.,8.,10.,20.])
    factors=(k*(b.H0/100)/Hm)[None,:]**2*beta[:,None]
    av=1/21
    result=[]
    for a0 in (.001,.005,.01):
        def rhs(u,v):
            a=np.exp(u)
            D,V=v.reshape((2,4,len(k)))
            return np.stack((V,a-.5*V-factors*D/a)).ravel()
        vals=[]
        for e in (1e-8,3e-10):
            rr=solve_ivp(rhs,np.log([a0,av]),np.zeros(2*4*len(k)),rtol=e,atol=1e-3*e,
                         method="DOP853")
            assert rr.success
            d=rr.y[:,-1].reshape(2,4,len(k))[0]
            d0=2/3*(av-3*a0+2*a0**1.5/np.sqrt(av))
            assert np.max(abs(d[-1]/d0-1.))<1e-8
            D=d[:3]/d0
            weighted=(1-f)*D[1]+f*D[2]
            vals.append((D,weighted))
        err=np.max(abs(vals[0][0]-vals[1][0]))
        assert err<3e-7
        D,S=vals[-1]
        result.append(dict(a_start=a0,z_final=20,k_h_per_Mpc=k.tolist(),joined_prompt_recoil=D.tolist(),
            weighted_separate=S.tolist(),joined_minus_weighted=(D[0]-S).tolist(),
            max_tolerance_control_defect=float(err)))
    return dict(input_sha256=[hashlib.sha256(f.read_bytes()).hexdigest() for f in (eos,*files)],
        recoil_number_fraction=float(f),q2_unscaled=q2.tolist(),results=result,
        scope="Newtonian EdS fixed-potential toy with zero initial perturbations, no shear, each supplied distribution closed with its own mean nonrelativistic adiabatic pressure coefficient. Nonlinear grouping sensitivity of a moment closure, not CLASS or cosmological transfer.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","joined","prompt","recoil"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,[a.joined,a.prompt,a.recoil]),indent=2,allow_nan=False))
