"""An additive fixed-background displacement control, not linear matter transfer."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import mp


def run(eos,files):
    b=Background(eos)
    ai=float(b.scale(.001))
    K=np.array([.1,.3,1.,2.,4.,8.,10.,20.])
    T=[np.loadtxt(f).T for f in files]
    cases=[]
    for z in (0,5,20):
        af=1/(1+z)
        controls=[]
        for size in (1800,3500):
            lo=min(np.log(q[0]) for q,_ in T)
            hi=max(np.log(q[-1]) for q,_ in T)
            q=np.exp(np.linspace(lo,hi,size))
            V=np.concatenate([b.streams(x*b.T0,mp,ai,af,7) for x in np.array_split(q,8)])
            zz=np.sinc(V[:,None]*K*b.H0/(100*np.pi))
            W=np.array([np.maximum(0.,np.nan_to_num(CubicSpline(x,y,extrapolate=False)(q),nan=0.))*q*q for x,y in T])
            N=np.trapezoid(W,q,axis=1)
            val=np.trapezoid(W[:,:,None]*zz[None,:,:],q,axis=1)/N[:,None]
            expected=np.average(val[1:],axis=0,weights=N[1:])
            err=np.max(np.abs(val[0]-expected))
            assert err<4e-7
            assert np.min(val)>-.218 and np.max(val)<=1.+1e-12
            controls.append((val,expected,N,err))
        drift=np.max(abs(controls[0][0]-controls[1][0]))
        assert drift<3e-5
        v,x,n,err=controls[-1]
        cases.append(dict(z=z,k_h_per_Mpc=K.tolist(),
              joined_prompt_recoil_displacement_characteristics=v.tolist(),
              weighted_two_component_characteristic=x.tolist(),
              recoil_number_fraction=float(n[2]/(n[1]+n[2])),
              max_absolute_joined_minus_split=float(err),max_refinement_drift=float(drift)))
    return dict(input_sha256=[hashlib.sha256(x.read_bytes()).hexdigest() for x in (eos,*files)],
         scope="Linear isotropic geodesic displacement characteristic (number-weighted j0 of k times free-streaming distance), without gravitational source terms or evolving perturbations. Exact linearity in the input number measure. Not a physical CLASS or total-matter transfer calculation.",cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    p.add_argument("joined",type=Path)
    p.add_argument("prompt",type=Path)
    p.add_argument("recoil",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,[a.joined,a.prompt,a.recoil]),indent=2,allow_nan=False))
