"""Use the actual ballistic parent occupation to bound the local pair gain."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m
from warm_pair_selective_loss_bound_reference import sha
from warm_scalar_pair_late_support_bound_reference import coefficient


def radial_bound(p,theta,n):
    # Exact p dp = (m+K) dK. The fitted entropy variable, and hence
    # the monotone FD occupation, have their extrema at the endpoints.
    K=p*p/(np.hypot(m,p)+m)
    a=np.linspace(0.,1.,n+1)
    f=expit(-theta[0,:,None]-theta[1,:,None]*a)
    z=K[:-1,None]+np.diff(K)[:,None]*a
    weights=np.diff(z,axis=1)*(2*m+z[:,:-1]+z[:,1:])/2.
    low=float(np.sum(np.minimum(f[:,:-1],f[:,1:])*weights))
    high=float(np.sum(np.maximum(f[:,:-1],f[:,1:])*weights))
    r,w=roots_legendre(44)
    X=K[:-1,None]+np.diff(K)[:,None]*(1.+r)/2.
    check=float(np.sum((m+X)*w*np.diff(K)[:,None]*expit(-theta[0,:,None]-theta[1,:,None]*(1+r)/2.)/2.))
    assert low<=check<=high
    return [low,check,high]


def run(shapes):
    out=[]
    for path in shapes:
        S=json.loads(path.read_text())["snapshots"][-1]
        T=S["T_visible_GeV"]
        q=np.asarray(S["lab_momentum_edges_GeV"])
        a,b=[np.asarray(S[k]) for k in ("entropy_variable_scalar","entropy_variable_fermion")]
        P=q[-1]
        Emin=m*m/(np.hypot(m,P)+P)
        Smax=float(expit(-np.min(b[0]+np.minimum(0.,b[1]))))
        z=np.min(a[0]+np.minimum(0.,a[1]))
        assert z>0.
        F0=float(np.exp(-z)/(-np.expm1(-z)))
        C=coefficient(P)
        H=np.sqrt(np.pi**2/45)/MPL*(T*.98)**2
        I=[radial_bound(q,b,n) for n in (4,16,64)]
        assert all(x[0]<=y[0] and x[-1]>=y[-1] for x,y in zip(I,I[1:]))
        # The p1 convolution has angular average 1/max(p1,k) <=
        # 1/p1. The p2 delta integral uses the same radial I/m.
        dose=32*C*P*P*4*I[-1][-1]**2/(128*Emin*(2*np.pi)**3*m*m*H)
        base=32*C*P**6*Smax*Smax/(128*Emin*(2*np.pi)**3*m*m*H)
        v=1./(1.+F0)
        assert 0.<dose<=base and dose<v
        out.append(dict(anchor_T=T,finite_incoming_pmax_GeV=P,source_shape_sha256=sha(path),
                        radial_lower_direct_upper_by_kinetic_panels_4_16_64=I,uniform_peak_based_dose=base,
                        radial_shape_dose_upper=dose,radial_to_uniform_dose=dose/base,initial_resolved_Bose_peak=F0,
                        global_Bose_ceiling=1./(v-dose)-1.,uniform_comparison_closes=bool(base<v)))
    return dict(cases=out,photon_floor_ratio=.98,spin_and_label_multiplier=32,
        scope="Same loss-only isolated ballistic comparison as the uniform peak bound; no inverse/redistribution/new parents and no inherited exterior scalar occupation. All produced scalar modes are allowed. The radial integral uses per-mode stable occupation rather than a number-density measure.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("shapes",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.shapes),indent=2,allow_nan=False))
