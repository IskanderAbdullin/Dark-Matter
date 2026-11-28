"""Apply the four existing scalar stress diagnostics to the new hard4 angular trajectory."""
import argparse
import json

import numpy as np
from scipy.optimize import brentq

from early_pole_daughter_bare_scalar_stress_control import run as bare, state
from early_pole_daughter_bare_scalar_frozen_source import run as frozen_source
from early_pole_daughter_retained_kinetic_stress import build, diagnostic, r
from early_pole_daughter_scalar_adiabatic_derivative import run as adiabatic
from early_pole_daughter_scalar_affine_shear import run as shear, check, paired_response


def spectator_control(C, arms):
    q=np.r_[C.q,C.spq]
    M=len(C.q)
    raw=[np.asarray(x['cold_mode_count_history'][-1]) for x in arms]
    def measures(zeta, verify=False):
        pair=[np.r_[n,zeta*C.spect] for n in raw]
        vals=[]
        for n in pair:
            x=state(q,n,C.lam,0.)
            E=np.hypot(q,x['effective_gap'])
            weight=q*q/E*(1-(q/E)**2/5)
            x['affine_shear']=float(n@weight)
            x['retained_gap_moment_fraction']=float((n[:M]@(1/E[:M]))/(n@(1/E)))
            if verify:
                independent=check(q,n,C.lam,0.)
                if abs(independent['affine_shear']/x['affine_shear']-1)>2e-14:
                    raise ValueError("Independent affine-shear replay mismatch.")
            if abs(x['pressure_to_energy']-1/3)>2e-13:
                raise ValueError("Spectator rescaling broke the conformal trace check.")
            vals.append(x)
        off=paired_response(q,*pair,*vals,M)
        return dict(factor=zeta,native=vals[0],local=vals[1],paired=off)
    root=brentq(lambda z:measures(z)['paired']['native_retained_response_fraction']-.5,1e-5,.1)
    tests=[measures(x,True) for x in (0.,1e-4,.001,root,.01,.1,1.,10.)]
    shares=[v['paired']['native_retained_response_fraction'] for v in tests]
    if any(a<=b for a,b in zip(shares[:-1],shares[1:],strict=True)) or abs(shares[3]-.5)>3e-11:
        raise ValueError("Spectator-share midpoint or sample ordering failed.")
    return dict(scope="Artificial spectator abundance rescaling at fixed terminal retained counts, bare mass zero.",
                equal_shear_spectator_factor=root,sweep=tests)


def run(archive, native, local, mesh):
    a=[r(native),r(local)]
    m=r(mesh)
    C,depth=build(archive,m)
    if depth!=4 or len(C.q)!=78:
        raise ValueError("Use the held hard4 angular78 pair.")
    v=diagnostic(C,depth,a)
    w=bare(archive,a,m)
    z=adiabatic(archive,a,w,m)
    p=shear(archive,a,z,m)
    return dict(scope="Fourth paired held scalar corner; no physical propagation or transfer result.",
                kinetic=v,bare=w,adiabatic=z,shear=p,
                frozen_initial_source=frozen_source(archive,m),spectator_stress_control=spectator_control(C,a))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','native','local','mesh'):
        p.add_argument(name)
    args=p.parse_args()
    print(json.dumps(run(**vars(args)),indent=2,allow_nan=False))
