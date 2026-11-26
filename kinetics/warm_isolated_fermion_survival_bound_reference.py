"""Conditional late survival of all dark fermions, not of one species alone."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from three_species_decay_onsager_reference import m,M,y
from warm_isolated_recycling_energy_envelope_reference import moments


def read(p):
    return json.loads(Path(p).read_bytes())


def one(T,E,N,omega,I,goal=.120):
    # For chi_i chibar_i -> two identical scalars, Cauchy after the direct
    # spin trace bounds sigma v_CM, then vM_lab/vM_CM<=1. Both charges
    # contribute: 2 n_plus n_minus <= n_total^2/2.
    KL,KH=y**4/(8*np.pi*np.array([M,m])**2)
    halfclock=.5*N*I
    h=min(E/(M*N),1.)
    assert h>0.
    EperN=E/N

    def g(z):
        high=min(z,max(0.,(EperN-m*z)/(M-m)))
        return max(KL*z*z,KL*(z-high)**2+KH*high*high)

    def exposure(z):
        assert 0<z<=1.
        if z<=h:
            base=(1/z-1/h)/KH
            z=h
        else:
            base=0.
        return base+quad(lambda t:1/g(t),z,1.,epsabs=1e-4,epsrel=3e-9,limit=300)[0]

    def lower(F,energy_weighted):
        d=(1+F)**2*halfclock
        coarse=1/(1+d*KH)
        if not energy_weighted:
            return coarse
        if d<exposure(h):
            return brentq(lambda z: exposure(z)-d,h,1.,xtol=1e-14)
        return 1/(1/h+KH*(d-exposure(h)))

    goalfrac=goal/omega
    critical_coarse=np.sqrt((1/goalfrac-1)/(KH*halfclock))-1
    critical_energy=np.sqrt(exposure(goalfrac)/halfclock)-1
    assert 0<critical_coarse<=critical_energy
    out=[]
    for F in (0.,.01,.1,1.,2.,2.5,3.):
        A,B=(lower(F,k) for k in (False,True))
        assert 0<A<=B*(1+1e-12)<=1.+1e-12
        out.append(dict(F_upper=F,coarse_fermion_rest_Omega_lower=A*omega,
                        energy_weighted_fermion_rest_Omega_lower=B*omega))
    return dict(anchor_T_GeV=T,core_start_fermion_number_GeV3=N,
                known_energy_upper_GeV4=E,conservative_initial_fermion_Omega=omega,
                integration_time_energy_assisted=I,coarse_minimum_cap_to_reach_goal=critical_coarse,
                energy_weighted_minimum_cap_to_reach_goal=critical_energy,controls=out)


def run(isolated,exterior,support,snapshots):
    A,B,C=map(read,(isolated,exterior,support))
    V=C["matched"][-1]["omegas_16q16_16q32_32q16_32q32"]
    rows=[]
    for i,path in enumerate(snapshots):
        S=read(path)["snapshots"][-1]
        e=B["cases"][i]
        a=A["cases"][i]
        T=S["T_visible_GeV"]
        assert abs(T/a["T_GeV"]-1.)<1e-12
        N,R=moments(S).T
        nr=N[-1]
        Er=R.sum()*(1+e["bounds"][-1]["added_energy_over_core"])
        # Assign both old and new counts to the heaviest stable species
        # to get a conservative lower *represented* stable rest proxy.
        extras=N.sum()*(e["known_new_count_over_core"]+e["known_old_count_over_core"])
        omega=V[i][1]/(1+extras/nr)
        I=a["bounds"][1]["clock"]["upper_integral_times_Hph"]/a["bounds"][1]["clock"]["photon_anchor_GeV"]
        rows.append(dict(q=e["q"],**one(T,Er,nr,omega,I)))
    return dict(scope="After the 10.442 MeV warm anchor, the specified isolated network "
                      "with arbitrary internal dark reentry/redistribution and arbitrary reverse "
                      "gains, but a uniform outgoing scalar occupation cap. No visible or "
                      "new fermion-number loss channel, no inference of the scalar cap.",
                target=.120,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("isolated","exterior","support"):
        p.add_argument(name)
    p.add_argument("snapshots",nargs=2)
    a=p.parse_args()
    print(json.dumps(run(a.isolated,a.exterior,a.support,a.snapshots),indent=2,allow_nan=False))
