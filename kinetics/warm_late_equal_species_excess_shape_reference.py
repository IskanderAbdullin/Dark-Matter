"""Match each daughter's convex energy separately while changing its packet shape."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from joint_prompt_two_moment_evolution_reference import kin
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import m, mp
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import excess, split_at_cut
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom

KEYS=("entropy_variable_scalar","entropy_variable_fermion")


def moments(D):
    x=get_mom(np.asarray(D["lab_momentum_edges_GeV"]),[D[k] for k in KEYS])
    return np.stack([x[0].sum(axis=1),x[2].sum(axis=1)])


def combine(Ds):
    ps=[np.asarray(x["lab_momentum_edges_GeV"]) for x in Ds]
    p=np.unique(np.r_[ps[0], ps[1]])
    assert ps[0][0]==ps[1][0] and ps[0][-1]==ps[1][-1]
    new=dict(Ds[0],lab_momentum_edges_GeV=p.tolist())
    for mass,k,q,D in zip((mp,m),KEYS,ps,Ds,strict=True):
        i=np.searchsorted(q,p[:-1],side="right")-1
        assert np.all(p[:-1]>=q[i]) and np.all(p[1:]<=q[i+1])
        a,b=np.asarray(D[k])
        en,src=kin(p,mass),kin(q,mass)
        d=np.diff(src)[i]
        new[k]=[(a[i]+b[i]*(en[:-1]-src[i])/d).tolist(),(b[i]*np.diff(en)/d).tolist()]
    for j,D in enumerate(Ds):
        np.testing.assert_allclose(moments(new)[:,j],moments(D)[:,j],rtol=3e-14)
    return new


def run(eos,archive,equal,factor=1.,extremes=False):
    saved=read(equal)
    C=saved["cut_GeV"]
    base=saved["cases"][0]["snapshot"]
    assert digest(base)==saved["cases"][0]["snapshot_sha256"]
    assert .5<factor<2.
    target=factor*np.asarray(excess(base,C))
    z=read(archive)
    S=z["snapshot"]
    assert digest(S)==z["snapshot_sha256"]==saved["source_sha256"]
    T=TabulatedRadiation(eos,S["T_visible_GeV"])

    def seed(p,w):
        D,*_=reactivate(S,p,w,cap=1.05)
        return split_at_cut(D,C)

    result=[]
    widths_to_try=((.0025,.0025),(.00025,.00025),(.00025,.0025),(.0025,.00025),(.001,.001))
    for widths in (widths_to_try[:2] if extremes else widths_to_try):
        centers=[]
        Ds=[]
        for j,w in enumerate(widths):
            f=lambda p: excess(seed(p,w),C)[j]/target[j]-1.
            a,b=.708,.717
            assert f(a)<0<f(b)
            center=brentq(f,a,b,xtol=1e-15)
            centers.append(center)
            Ds.append(seed(center,w))
        D=combine(Ds)
        U=moments(D)
        rho=U[1].sum()
        H=np.sqrt((T.state(1.)[1]+rho)/3.)/MPL
        parts=np.asarray(excess(D,C))
        dw=float(max(abs(parts/target-1.)))
        dm=float(np.max(abs(U/moments(base)-1.)))
        assert dw<3e-11 and dm<4e-14
        rates=[]
        for order in (8,12):
            g=integrate(D,C,order,True)/(rho*H)
            rates.append(dict(order=order,signed=float(g[3]),gross=float(g[1]+g[2])))
        assert abs(rates[0]["signed"]-rates[1]["signed"])<2e-5*rates[-1]["gross"]
        result.append(dict(widths_GeV=widths,centers_GeV=centers,W_match_error=dw,moment_match_error=dm,
                           W_over_rho= (parts/rho).tolist(),rates=rates,hash=digest(D),snapshot=D))
    return dict(cut_GeV=C,archive_hash=z["snapshot_sha256"],reference_hash=digest(base),
                target_factor=factor,target_W_over_rho=(target/moments(base)[1].sum()).tolist(),
                scope="Frozen same-support statistical prompt controls; independent species packet shapes "
                "exactly restrict to a common grid before angular integration.",cases=result)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","equal"):
        p.add_argument(k,type=Path)
    p.add_argument('--target-factor',type=float,default=1.)
    p.add_argument('--extremes-only',action='store_true')
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive,a.equal,a.target_factor,a.extremes_only),indent=2,allow_nan=False))
