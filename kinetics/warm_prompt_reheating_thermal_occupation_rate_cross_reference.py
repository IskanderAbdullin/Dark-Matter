"""Pointwise thermal identity and quantum-rate controls on the dilute expanding arm."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from archived_scalar_pair_reaction_reference import gross,points
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_heavy_gain_probe_reference import rates
from warm_prompt_reheating_signed_slow_closure_reference import pair_mb


def run(path):
    data=json.loads(path.read_text())
    if "arms" in data or "controls" in data:
        if "arms" in data:
            r=data["arms"][0]["rows"]
            key,lab="T_dark","peak"
        else:
            r=data["controls"][-1]["rows"]
            key,lab="dark_T","phase_factor"
        take=(r[0],max(r,key=lambda q:max(q[lab])),r[-1])
        seq=[(j,dict(T_dark_GeV=x[key],max_phase_occupation=x[lab])) for j,x in enumerate(take)]
    else:
        seq=[(j,data["paths"][1]["rows"][j]) for j in (0,4,10)]
    out=[]
    for j,p in seq:
        T=p["T_dark_GeV"]
        pk=np.array(p["max_phase_occupation"])
        eta=np.log(pk)+np.array((m,mp,M))/T
        assert abs(eta[2]-eta[0]-eta[1])<2e-12
        theta=np.exp(eta)
        def f(e):
            x=e/T-eta[1]
            return np.exp(-x)/(-np.expm1(-x))
        g=lambda e:expit(eta[0]-e/T)
        lo=np.array([1/((1+pk[0])*(1+pk[2]))**2,1/(1+pk[2])**2,1/(1+pk[0])**2])
        hi=np.array([1.,1/(1-pk[1])**2,1/(1-pk[1])**2])
        C=pair_mb(T,64)*np.array([theta[0]**2,theta[1]**2,theta[0]**2])
        tests=[]
        orders=((50,3,15,10),(70,4,21,16)) if max(pk)<.4 else ((62,3,23,18),(85,4,34,25))
        for n,rad,ang,az in orders:
            edge=np.r_[0.,np.geomspace(.0004,max(70*T,2.),n)]
            pp,_,u=points(edge,m,g,rad)
            ss,_,v=points(edge,mp,f,rad)
            x=np.array((rates(pp,4*u,m,ang).sum(),rates(ss,v,mp,ang).sum()))
            _,z,_=gross(edge,m,mp,g,f,False,rad,ang,az,10,False,exchange_mass=M)
            _,w,_=gross(edge,mp,m,f,g,True,rad,ang,az,10,False,exchange_mass=M)
            d=w[1]/z[1]/np.exp(2*(eta[1]-eta[0]))-1
            _,mb,_=gross(edge,m,mp,lambda e:np.exp(eta[0]-e/T),np.zeros_like,
                           False,rad,ang,az,10,False,exchange_mass=M)
            # The bound is pointwise; use the identical positive nodes for
            # the ratio so an ultra-dilute correction is not masked by the
            # separate continuum MB normalization error.
            h=z[1]/mb[0]
            assert lo[2]-1e-13<=h<=hi[2]+1e-13 and abs(d)<.003,(j,T,h,float(d),lo[2],hi[2])
            tests.append(dict(radial=rad,angle=ang,edge_count=n,
                              light_to_scalar_forward_multiplier_over_MB=float(h),
                              matched_lab_MB_over_invariant_MB=float(mb[0]/C[2]),
                              light_to_scalar_reverse_forward_ratio_relative_identity_defect=float(d),
                              heavy_forward_without_final_pauli_over_MB=(x/C[:2]).tolist(),
                              conditional_forward_multiplier_bounds=[lo.tolist(),hi.tolist()]))
        out.append(dict(sample=j,T_GeV=T,controls=tests))
    return dict(scope="Pointwise thermal rate check, not time integration. Exact FD/BE incoming and final factors on the light/scalar pair; "
                      "heavy forward test omits heavy final Pauli. Sampling convergence tested independently.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("sensitivity",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.sensitivity),indent=2,allow_nan=False))
