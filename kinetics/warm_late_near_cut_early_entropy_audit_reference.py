"""Read the first rejected small-step root; do not accept it."""
import json
import math
import argparse
from pathlib import Path

import numpy as np

from joint_prompt_two_moment_evolution_reference import Collision, rule
from joint_resonant_centered_stiff_reference import advance
from joint_resonant_frw_split_reference import pack
from warm_compact_absolute_prompt_reference import LocalAbsolutePrompt
from warm_corrected_late_prompt_cut_aligned_reference import predictor
from warm_light_occupation_stage_audit_reference import read


def entropy_line(op, u, v, order):
    """Integrate ds = theta * f'(theta) dtheta without large entropy subtraction."""
    a,w=rule(order)
    t=(a+1)/2
    records=[]
    for c,ac,bc in zip(op.k.c,op.arrays(u),op.arrays(v),strict=True):
        dtype=np.longdouble
        wq=c.vol.astype(dtype)[:,None]*c.w.astype(dtype)
        uc=c.u.astype(dtype)
        da,db=(bc.astype(dtype)-ac.astype(dtype))
        z0=ac[0].astype(dtype)[:,None]+ac[1].astype(dtype)[:,None]*uc
        dz=da[:,None]+db[:,None]*uc
        z=z0[:,:,None]+t*dz[:,:,None]
        positive=np.exp(-abs(z))
        if c.g==1:
            assert np.all(z>0.)
            variance=positive/(-np.expm1(-z))**2
        else:
            variance=positive/(1+positive)**2
        out=-(wq*dz)*(z*variance @ (w/2))
        records.append(out)
    x=np.concatenate([r.ravel() for r in records])
    return math.fsum(x), math.fsum(abs(x)), [math.fsum(r.ravel()) for r in records]


def stable_total(op, theta):
    # Two numerically complementary arguments to log(1-exp(-z)).
    # Retain long-double arithmetic until after differencing.
    out=np.longdouble(0)
    for c, ac in zip(op.k.c,op.arrays(theta),strict=True):
        z=ac[0,:,None].astype(np.longdouble)+ac[1,:,None].astype(np.longdouble)*c.u
        y=np.exp(-z)
        if c.g==1:
            f=y/(-np.expm1(-z))
            v=np.empty_like(z)
            soft=z<np.log(2)
            v[soft]=-np.log(-np.expm1(-z[soft]))
            v[~soft]=-np.log1p(-y[~soft])
        else:
            f=np.where(z>=0,y/(1+y),1/(1+np.exp(z)))
            v=np.logaddexp(0,-z)
        out+=np.sum(c.vol.astype(np.longdouble)[:,None]*c.w*(z*f+v))
    return out


def run(p):
    raw=read(p/"step-0.json")
    S=raw["snapshot"]
    D=read(p/"summary.json")
    assert D["plan"]==raw["plan"]
    H=D["H_GeV"]
    op=LocalAbsolutePrompt(Collision(np.asarray(S["lab_momentum_edges_GeV"])),H,[],True,stream_local=True)
    u=pack([S["entropy_variable_"+k] for k in ("scalar","fermion")])
    h=1e-11/6
    for j in range(6):
        rejected={}
        def inspect(v,s,s0):
            if s>=s0:
                return True
            U=op.moment(v)[0]
            U0=op.moment(u)[0]
            F=op.collision(v)[0]
            err=(U-U0)-h*F
            diffs=[entropy_line(op,u,v,n) for n in (8,16,32)]
            integrated=diffs[-1][0]
            stable=float(stable_total(op,v)-stable_total(op,u))
            lower=math.fsum(v*(U-U0))
            assert integrated>0. and stable>0. and lower>0.
            assert max(abs(a[0]-integrated) for a in diffs)<1e-7*integrated
            assert abs(stable-integrated)<4e-5*integrated
            rejected.update(index=j+1,h=h,delta_entropy_float=float(s-s0),
                total_entropy=float(s0),ulp=float(np.spacing(s0)),line_controls=diffs,
                independent_stabilized_absolute_difference=stable,
                terminal_collisional_entropy_times_dt=h*math.fsum(v*F),
                approximate_root_entropy_residual=math.fsum(v*err),
                thermodynamic_concavity_lower=lower,anchor_hash=raw["hash"])
            return False
        try:
            u,_=advance(op,u,h,S["T_visible_GeV"],initial_theta=predictor(op,u,h),
                        vacuum_scaled=True,balanced_merit=True,entropy_guard=inspect)
        except AssertionError:
            if not rejected:
                raise
            return rejected
    raise RuntimeError("Did not reproduce the initial rejection")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("original_phase_folder",type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.original_phase_folder),indent=2,allow_nan=False))
