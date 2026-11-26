"""A conditional one-sided shell regularity criterion for a peak cap."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from three_species_decay_onsager_reference import m
from warm_light_occupation_stage_audit_reference import digest, read
from warm_shell_occupation_feasibility_reference import independent_moments, kinetic, volume


def minimum_number(a, b, c, slope):
    """Sharp radial number lower bound for peak >= c and |df/dp| <= slope."""
    d = min(b-a, c/slope) if slope>0 else b-a
    t = slope*d/c
    return c*d*(a*a*(1-t/2)+a*d*(1-2*t/3)+d*d*(1/3-t/4))


def threshold(a, b, c, n):
    if not 0<n<c*volume(a,b):
        raise ValueError("need a positive number below the flat-cap average")
    # Solve in log slope for very dilute shells, where the critical slope
    # may be many decades above the scale c/(b-a).
    s0=c/(b-a)
    def f(log_slope):
        return math.log(minimum_number(a,b,c,math.exp(log_slope))) - math.log(n)
    low=math.log(s0)
    while f(low)<=0:
        low-=4.
    high=math.log(s0)
    while f(high)>=0:
        high+=8.
    return math.exp(brentq(f,low,high,xtol=2e-13))


def verify():
    from scipy.integrate import quad
    above=[]
    test=[(0.,1.,.5,1.),(.1,.7,.5,1.),(.1,.7,.5,.1),(.1,.7,.5,0.)]
    for a,b,c,L in test:
        base=minimum_number(a,b,c,L)
        for x in np.linspace(a,b,9):
            points=[x] if L==0. else [t for t in (x-c/L,x,x+c/L) if a<t<b]
            exact=quad(lambda y: y*y*max(c-L*abs(y-x),0.),a,b,points=points)[0]
            assert exact>=base*(1-3e-15)
            if x==a:
                assert abs(exact/base-1)<3e-15
        above.append(base)
    assert abs(above[0]-(.5**4)/12)<5e-18
    return above


def run(folder):
    checks=verify()
    rows=[]
    for q in (16,32):
        frames=read(folder/f'warm-light-occupation-next-q{q}-tape-2026-11-26.json.gz')
        frames+=[z["checkpoint"] for z in
                  read(folder/f'warm-observed-suboperator-q{q}-tape-2026-11-26.json.gz')["blocks"]]
        stages=[]
        for x in frames:
            s=x["snapshot"]
            assert x["snapshot_sha256"]==digest(s)
            P=np.asarray(s["lab_momentum_edges_GeV"])
            uv=np.asarray(s["entropy_variable_fermion"])
            caps=[]
            for j,(a,b,(u,v)) in enumerate(zip(P[:-1],P[1:],uv.T,strict=True)):
                n,_=independent_moments(a,b,m,[u,v])
                if n<=0:
                    continue
                L=threshold(a,b,.5,n)
                reconstructed_bound=float(abs(v)*expit(-min(u,u+v))*(b/np.hypot(b,m))/(kinetic(b,m)-kinetic(a,m)))
                assert abs(minimum_number(a,b,.5,L)/n-1)<5e-12
                assert reconstructed_bound<L
                caps.append(dict(shell_index=j, slope_GeV_inverse=L,
                                 slope_in_normalized_linear_shell=L*(b-a),
                                 reconstructed_slope_GeV_inverse_upper=reconstructed_bound))
            stages.append(dict(T_GeV=s["T_visible_GeV"], snapshot_sha256=x["snapshot_sha256"],
                               minimum_dimensional_threshold=min(caps,key=lambda z:z["slope_GeV_inverse"]),
                               minimum_normalized_threshold=min(caps,key=lambda z:z["slope_in_normalized_linear_shell"]),
                               max_reconstructed_to_critical_slope=max(v["reconstructed_slope_GeV_inverse_upper"]/
                                                       v["slope_GeV_inverse"] for v in caps)))
        rows.append(dict(q=q,stages=stages))
    return dict(occupation_cap=.5, analytical_cone_controls=checks,cases=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tapes",type=Path,default=Path("results"))
    print(json.dumps(run(parser.parse_args().tapes),indent=2,allow_nan=False))
