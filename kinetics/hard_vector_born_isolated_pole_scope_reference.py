"""Toy free-quark Born weight in a disjoint high-mass interval around the pole."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from born_thermal_quark_hard_portal_reference import UEDGES, gauss, local, thermal_trace_over_vacuum
from conditional_relic_normalization_reference import conversion
from disjoint_hard_thermal_measured_source_reference import WINDOWS
from measured_partial_hadron_reheat_fold_reference import Fold
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import mp
from timelike_portal_onsager_coverage_reference import V, width

RHO=V*width
BANDS=np.array([8.5,V-.05,V+.05,11.5])
PANELS=np.unique(np.clip(np.r_[BANDS,V+width*np.array([-200,-60,-20,-8,-3,-1,0,1,3,8,20,60,200])],BANDS[0],BANDS[-1]))


def at(T,a,b,coordinate,n):
    u,du=gauss(UEDGES,n)
    def z(E,dEdv):
        x=(E/T)[:,None]
        e=x+u
        w=4*dEdv[:,None]/T*du*np.sqrt(u*(2*x+u))/np.expm1(e)*thermal_trace_over_vacuum(x,u)
        return local(T,(x[:,0],w.sum(axis=1)[None,:],(w*e).sum(axis=1)[None,:]))[0]
    edges=np.unique(np.clip(PANELS,a,b))
    if coordinate=="arctan-s":
        v,w=gauss(np.arctan((edges**2-V**2)/RHO),n)
        E=np.sqrt(V*V+RHO*np.tan(v))
        return z(E,w*RHO/(2*E*np.cos(v)**2))
    if coordinate=="log-M":
        v,w=gauss(np.log(edges),n)
        E=np.exp(v)
        return z(E,w*E)
    raise ValueError(coordinate)


def integrate(F,fac,Tmin,Tmax,n,coordinate):
    aa,bb=np.log([Tmin,Tmax])
    edge=np.unique(np.clip(np.r_[aa,bb,F.edges],aa,bb))
    t,w=gauss(edge,n)
    out=np.zeros((3,2))
    for v,dv in zip(t,w):
        T=np.exp(v)
        s,r=F.densities(v)
        H=np.sqrt(r/3)/MPL
        dt=dv*(3+F.bg.sg.derivative()(v))/(3*H)
        for j,(a,b) in enumerate(zip(BANDS[:-1],BANDS[1:])):
            C=at(T,a,b,coordinate,n)
            out[j,0]+=fac*dt/s*(C[0,0]+2*C[1,0])
            out[j,1]+=dt*(F.sf/s)**(4/3)/F.rf*C[:,1].sum()
    return out


def run(eos):
    # Keep this mass band disjoint from both measured windows and the
    # kinematic support of the previous sub-pole (x≤24, T≤.345) diagnostic.
    assert max(WINDOWS.ravel())<24*.345<BANDS[0]
    F=Fold(eos,(8,21,16),6,upper=.50)
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    rows=[]
    for T0,T1 in ((.25,.32236),(.30,.33577)):
        ref=integrate(F,f,T0,T1,12,"log-M")
        theta=integrate(F,f,T0,T1,12,"arctan-s")
        low=integrate(F,f,T0,T1,8,"arctan-s")
        e=float(np.max(abs(theta/ref-1.)))
        d=float(np.max(abs(low/theta-1.)))
        assert np.all(theta>0.) and max(e,d)<5e-6
        tot=theta.sum(axis=0)
        rows.append(dict(T_limits_GeV=[T0,T1],isolated_band_partition_count_and_heat=theta.tolist(),
                         total_same_band=tot.tolist(),share_within_pole_window=(theta[1]/tot).tolist(),
                         fraction_of_reference_optimistic_count_target=float(tot[0]/.1188),
                         largest_coordinate_relative_difference=e,largest_order_relative_difference=d))
    return dict(mass_panels_GeV=BANDS.tolist(),vector_width_GeV=width,scope_max_previous_M_GeV=24*.345,
                EOS_sha256=hashlib.sha256(eos.read_bytes()).hexdigest(),
                columns=["Omega h^2 if all convert","dark energy floor at 1 MeV / radiation"],
                controls=rows,scope="Isolated high-mass massless three-flavor Born proxy with the same Breit-Wigner propagator and imposed radiation clock. No measured high-energy spectrum, heavy flavors, other disjoint intervals, coupled evolution, full QCD correction, or fitted reheating root.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos),indent=2,allow_nan=False))
