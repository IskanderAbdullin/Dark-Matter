"""Independent invariant-mass panels for the measured-windows exclusion."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from born_thermal_quark_hard_portal_reference import local
from partial_hard_vector_nlo_reheat_reference import Correction, Fold, WINDOWS, XEDGES, interp, table


def independent(T,values,mode,n,lower_x=None):
    X,U,C=values
    p,w=leggauss(n)
    ua,ub=U[:-1,None],U[1:,None]
    u=((ua+ub+(ub-ua)*p)/2).ravel()
    du=((ub-ua)*w/2).ravel()
    evaluator=interp(X,U,C,mode)(u)
    assert lower_x is None or X[0] <= lower_x < X[-1]
    xe=X if lower_x is None else np.unique(np.r_[lower_x,X[X>lower_x]])
    mx=T*xe
    panels=np.unique(np.clip(np.r_[mx,WINDOWS.ravel()],mx[0],mx[-1]))
    result=np.zeros((1+len(WINDOWS),2,2))
    for a,b in zip(panels[:-1],panels[1:]):
        mass=((a+b)+(b-a)*p)/2
        x=mass/T
        # Jacobian converts independent dM back to the factor dx in local.
        dx=(b-a)*w/(2*T)
        E=x[:,None]+u
        k=np.sqrt(u*(2*x[:,None]+u))
        f=4*dx[:,None]*du*k*evaluator(x)/np.expm1(E)
        center=(a+b)/2
        j=np.flatnonzero((WINDOWS[:,0]<center)&(center<WINDOWS[:,1]))
        assert len(j)<=1
        i=0 if not len(j) else int(j[0])+1
        result[i]+=local(T,(x,f.sum(axis=1)[None,:],(f*E).sum(axis=1)[None,:]))[0]
    return result


def run(eos, native, parent_x_max=12.):
    if parent_x_max not in (12.,13.,14.):
        raise ValueError("Supported Born upper edges are 12, 13, and 14.")
    parent=np.r_[XEDGES,np.arange(13.,parent_x_max+1.)]
    values=table(native,parent[-1])
    low,high=values[0][[0,-1]]
    first_cross=WINDOWS[0,0]/high
    last_cross=WINDOWS[0,1]/high
    second_cross=WINDOWS[1,0]/high
    assert first_cross<second_cross
    fold=Fold(eos,(6,16,12),4,upper=.5)
    tests=[]
    for mode in ("linear","pchip"):
        solver=Correction(fold,values,.25,8,mode,parent)
        probes=[first_cross-1e-6,first_cross+1e-6,.32517560695,.33730194806,.40,last_cross+1e-5,.47]
        if second_cross<.48:
            probes.extend([second_cross-1e-6,second_cross+1e-6])
        for T in sorted(set(probes)):
            proposed=solver.rates(T)
            V=independent(T,values,mode,16)
            V2=independent(T,values,mode,24)
            outside,regions=V2[0],V2[1:]
            blocked=regions.sum(axis=0)
            total=outside+blocked
            assert np.all(regions[1]==0.) if T<second_cross else np.all(regions[1]>0.)
            ref=np.array([total,outside])
            assert np.all(outside>0) and np.all(blocked>=0)
            assert np.all(blocked==0.) if T<first_cross else np.all(blocked>0.)
            err=np.max(abs((proposed-ref)/total))
            convergence=np.max(abs(V-V2)/total)
            assert err<1e-6 and convergence<2e-7,(T,mode,err,convergence)
            tests.append(dict(interpolation=mode,T_GeV=T,
                local_blocked_fraction_count_or_energy_by_stable_heavy=(blocked/total).tolist(),
                local_each_window_blocked_fraction_separately=(regions/total).tolist(),
                max_deviation_repartition_normalized_to_unsplit=float(err),
                independent_quadrature_relative_deviation=float(convergence)))
    return dict(mass_over_T_interval=[float(low),float(high)], first_BES_entry_T_GeV=first_cross,
        first_BES_exit_at_top_T_GeV=last_cross,second_BES_entry_T_GeV=second_cross,
        strip_sha256=hashlib.sha256(Path(native).read_bytes()).hexdigest(),
        controls=tests,scope="Checks invariant-mass windows and numerical partition only, not QCD validity or the conversion target.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("native")
    p.add_argument("--born-x-max",type=float,default=12.)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.native,a.born_x_max),indent=2,allow_nan=False))
