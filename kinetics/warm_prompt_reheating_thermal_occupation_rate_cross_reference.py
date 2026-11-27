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
    out=[]
    for j in (0,4,10):
        p=data["paths"][1]["rows"][j]
        T=p["T_dark_GeV"]
        pk=np.array(p["max_phase_occupation"])
        eta=np.log(pk)+np.array((m,mp,M))/T
        assert abs(eta[2]-eta[0]-eta[1])<2e-12
        theta=np.exp(eta)
        f=lambda e:1/np.expm1(e/T-eta[1])
        g=lambda e:expit(eta[0]-e/T)
        lo=np.array([1/((1+pk[0])*(1+pk[2]))**2,1/(1+pk[2])**2,1/(1+pk[0])**2])
        hi=np.array([1.,1/(1-pk[1])**2,1/(1-pk[1])**2])
        C=pair_mb(T,64)*np.array([theta[0]**2,theta[1]**2,theta[0]**2])
        tests=[]
        for n,rad,ang,az in ((50,3,15,10),(70,4,21,16)):
            edge=np.r_[0.,np.geomspace(.0004,max(70*T,2.),n)]
            pp,_,u=points(edge,m,g,rad)
            ss,_,v=points(edge,mp,f,rad)
            x=np.array((rates(pp,4*u,m,ang).sum(),rates(ss,v,mp,ang).sum()))
            _,z,_=gross(edge,m,mp,g,f,False,rad,ang,az,10,False,exchange_mass=M)
            _,w,_=gross(edge,mp,m,f,g,True,rad,ang,az,10,False,exchange_mass=M)
            d=w[1]/z[1]/np.exp(2*(eta[1]-eta[0]))-1
            h=z[1]/C[2]
            assert lo[2]<h<hi[2] and abs(d)<.003
            tests.append(dict(radial=rad,angle=ang,edge_count=n,
                              light_to_scalar_forward_multiplier_over_MB=float(h),
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
