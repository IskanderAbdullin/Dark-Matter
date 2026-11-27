"""Full BE/FD prompt-return heat gap on held thermal continuation states."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre,expit

from archived_collisionless_streaming_reference import Background
from factorized_prompt_return_graph_reference import C
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp,P,Es
from warm_prompt_reheating_expanding_slow_reference import visible


def geometry(T,n,j):
    u,w=roots_legendre(n)
    x,v=roots_legendre(j)
    z=(u+1)*np.sqrt(67.)/2
    E=M+T*z*z
    dE=T*z*np.sqrt(67.)*w
    p=np.sqrt((E-M)*(E+M))
    hi=(E*Es+p*P)/M
    lo=(Es*Es+(p*mp/M)**2)/hi
    delta=np.log(hi/lo)
    e=lo[:,None]*np.exp(delta[:,None]*(1+x)/2)
    return E,dE,e,delta[:,None]*e*v/2


def thermal(T,eta,n,j,tilt=None):
    E,dE,e,w=geometry(T,n,j)
    Tphi=T if tilt is None else T*(1+tilt)
    f=1/np.expm1(e/Tphi-eta[1])
    g=expit(eta[0]-(E[:,None]-e)/T)
    A=w*f*g
    B=w*(1+f)*(1-g)
    b=B.sum(axis=1)
    mu=(B*e).sum(axis=1)/b
    dx=e-mu[:,None]
    if tilt is None:
        r=expit(eta[0]+eta[1]-E/T)
        return C*np.sum(dE*r*np.sum(B*dx*dx,axis=1))/T**2
    # The exact signed scalar energy rate in the prompt projection.
    # Centering E_phi removes cancellation of O(1) offsets.
    return -C*np.sum(dE*b/(A.sum(axis=1)+b)*np.sum(A*dx,axis=1))


def capacity(T,eta,mass,g,bose,n):
    z,w=roots_legendre(n)
    bounds=np.r_[0.,np.geomspace(1e-5,77.,29)]
    q=T*((bounds[:-1,None]+bounds[1:,None])/2+np.diff(bounds)[:,None]*z/2)
    E=np.hypot(q,mass)
    f=1/np.expm1(E/T-eta) if bose else expit(eta-E/T)
    k=g*w*np.diff(bounds)[:,None]*T*q*q*f*(1+f if bose else 1-f)/(4*np.pi**2)
    av=(k*E).sum()/k.sum()
    return (k*(E-av)**2).sum()/T**2


def run(eos,history):
    bg=Background(eos)
    data=json.loads(history.read_text())
    if "arms" in data:
        samples=[dict(dark_T=r["T_dark"],visible_T=r["T_visible"],
                      phase_peak=r["peak"],invariants=r["invariant_yields"],
                      d_invariants_over_Hs=np.sum(r["dark_and_visible_count_sources"],axis=0).tolist())
                 for r in data["arms"][0]["rows"]]
    else:
        samples=data["history"][-1]["samples"]
    out=[]
    for i,p in enumerate(samples):
        T=p["dark_T"]
        eta=np.log(p["phase_peak"])+np.array((m,mp,M))/T
        rv,sv=visible(bg,p["visible_T"])
        H=np.sqrt(rv/3)/MPL
        levels=[]
        for n,j in ((70,42),(104,68)):
            K=thermal(T,eta,n,j)
            z=2e-5
            check=-(thermal(T,eta,n,j,z)-thermal(T,eta,n,j,-z))/(2*z*T)
            assert abs(check/K-1)<8e-7
            Cv=np.array([capacity(T,eta[0],m,4,False,n//4),capacity(T,eta[1],mp,1,True,n//4)])
            rate=K*np.sum(1/Cv)/H
            rate_slow=np.max(abs(np.array(p["d_invariants_over_Hs"])/np.array(p["invariants"])))
            levels.append(dict(orders=[n,j],conductance=K,finite_difference_relative_check=float(abs(check/K-1)),
                               heat_gap_over_visible_H=float(rate),largest_fractional_slow_number_rate_over_H=float(rate_slow),
                               heat_to_largest_slow_rate_ratio=float(rate/max(1.,rate_slow))))
        assert abs(levels[0]["heat_gap_over_visible_H"]/levels[1]["heat_gap_over_visible_H"]-1.)<6e-6
        out.append(dict(i=i,visible_T=p["visible_T"],dark_T=T,peak_occupation=p["phase_peak"],levels=levels))
    return dict(scope="Local linear directed energy-exchange gap of the prompt-return, fixed-count, two-light-species model. "
                      "Full FD and BE for the two outgoing legs. Heavy heat capacity and off-manifold distortions "
                      "are not included. Comparison to selected slower number rates is not a spectral gap "
                      "of the full collision operator.",rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("history",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.history),indent=2,allow_nan=False))
