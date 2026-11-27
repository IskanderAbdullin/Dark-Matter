"""Phase-resolved pressure and velocity moments for the local coupled endpoints."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def observable(s,order):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    B,_,F,_=snapshot_functions(s)
    z,w=roots_legendre(order)
    r=[]
    for mass,deg,fun in ((mp,1,B),(m,4,F)):
        poles=mass*np.array([.1,1.,10.])
        edges=np.unique(np.r_[p,poles[(poles>p[0])&(poles<p[-1])]])
        dp=np.diff(edges)[:,None]
        q=edges[:-1,None]+dp*(z+1)/2
        e=np.hypot(q,mass)
        v2=(q/e)**2
        phase=deg/(4*np.pi**2)*dp*w*q*q*fun(e)
        N=phase.sum()
        E=np.sum(phase*e)
        pressure=np.sum(phase*e*v2)/3
        trace=np.sum(phase*(mass*mass/e))
        assert abs((3*pressure+trace)/E-1.)<2e-13
        r.append([N,E,pressure,trace,np.sum(phase*v2)/N,
                  *[np.sum(phase[q<=t])/N for t in poles]])
    r=np.asarray(r)
    N,_,E=get_mom(p,[s["entropy_variable_"+k] for k in ("scalar","fermion")])
    np.testing.assert_allclose(r[:,:2],np.stack((N.sum(axis=1),E.sum(axis=1)),axis=1),rtol=2e-11,atol=0.)
    return r


def run(anchor,*cases):
    r=read(anchor)
    s=r["snapshot"]
    assert digest(s)==r["hash"]
    initial=observable(s,48)
    base=None
    out=[]
    for p in cases:
        a=read(p)
        s=a["daughter_snapshot"]
        row=a["row"]
        assert digest(s)==row["daughter_hash"] and row["initial_daughter_hash"]==r["hash"]
        x,y=[observable(s,n) for n in (32,48)]
        np.testing.assert_allclose(x,y,rtol=2e-11,atol=0.)
        if base is None:
            base=y
        out.append(dict(case=p.parent.name,hash=row["daughter_hash"],a=row["c"],quantities=y.tolist(),
                        fractional_shift_vs_first=np.divide(y-base,base,out=np.zeros_like(y),where=base>0.).tolist(),
                        cold_number_fraction_absolute_shift_vs_first=(y[:,5:]-base[:,5:]).tolist()))
    return dict(scope="Instantaneous kinematic phase-space moments only; no post-window collisionless history, thermal-equivalent temperature, or CLASS prediction.",
                initial=initial.tolist(),initial_hash=r["hash"],species=["scalar","light fermion"],
                columns=["N","rho","pressure","rho_minus_3P","number_weighted_v_squared","number_fraction_p_below_.1m",
                         "number_fraction_p_below_m","number_fraction_p_below_10m"],rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("anchor",type=Path)
    p.add_argument("cases",type=Path,nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.anchor,*a.cases),indent=2,allow_nan=False))
