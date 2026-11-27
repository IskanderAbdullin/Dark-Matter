"""Two-stage L-stable collision step with stage-local parent elimination."""
import sys

import numpy as np

from joint_resonant_centered_stiff_reference import advance


GAMMA=1.-1./np.sqrt(2.)


def stability_test():
    """Check the tableau on a scalar relaxation, without a positivity claim."""
    g=GAMMA
    np.testing.assert_allclose([g+(1.-g), (1.-g)*g+g], [1.,.5], rtol=3e-15)
    for z in (-.01,-2.,-50.,-1e8):
        y1=1./(1.-g*z)
        y2=(1.+(1.-g)*z*y1)/(1.-g*z)
        expected=(1.+(1.-2.*g)*z)/(1.-g*z)**2
        np.testing.assert_allclose(y2,expected,rtol=2e-14)
    # The late relaxation remainder is not universally positivity preserving.
    assert expected<0. and abs(expected)<1e-7


def step(op,u,h,T,guide,fin,log,tol):
    g=GAMMA
    U0,S0,_=op.moment(u)
    h0=op.old.copy() if fin else None
    old_parent_entropy=op.heavy_entropy(h0) if fin else 0.
    def guard(v,after,before):
        if fin:
            z=op.occupation(v)
            return (np.all((0.<=z)&(z<=1.))
                    and after+op.heavy_entropy(z)>=before+old_parent_entropy-8e-15*(before+old_parent_entropy))
        return after>=before-8e-15*before
    options=dict(vacuum_scaled=True,balanced_merit=True,
                 entropy_guard=guard,residual_tolerance=tol)
    if fin:
        op.dt=g*h
    ctx={}
    q1,r1=advance(op,u,g*h,T,initial_theta=u+g*(guide-u),observer=log,
                  centered_diagnostic=ctx,**options)
    F=op.centered(ctx["centered"],ctx["reference"],ctx["parent_ratio"])[0]
    h1=None
    affine_bounds=None
    if fin:
        h1=op.occupation(q1)
        _,_,_,_,A,B=op.reactions(q1)
        strength=g*h*op.relax_factor
        inc=strength*(A-(A+B)*h0)/(1.+strength*(A+B))
        # The affine RK target is algebra, not a physical occupation.
        # Only the actual two implicit stages must lie in [0,1].
        op.old=h0+(1.-g)/g*inc
        affine_bounds=[float(op.old.min()),float(op.old.max())]
    rhs=U0+(1.-g)*h*F
    a=rhs.reshape(-1,2)
    # This is an affine RK target, not an intermediate physical
    # daughter reconstruction. Leave the actual target unchanged.
    nonreal=int(np.sum((a[:,0]<0.)|(a[:,1]<0.)|(a[:,1]>a[:,0])))
    print("sdirk2 stage-2 affine target outside the simple cone in",
          nonreal,"cells",file=sys.stderr,flush=True)
    rejected=[]
    for label,guess in (("stage_shift",q1+(1.-g)*(guide-u)),("stage_one",q1),("guide",guide)):
        try:
            q2,r2=advance(op,u,g*h,T,initial_theta=guess,target_moments=rhs,observer=log,**options)
            break
        except (RuntimeError,ValueError,AssertionError,FloatingPointError,RuntimeWarning) as e:
            rejected.append((label,str(e) or type(e).__name__))
            print("sdirk2 stage-2 seed rejected:",rejected[-1],file=sys.stderr,flush=True)
    else:
        raise RuntimeError(f"no admissible SDIRK2 stage-2 solution: {rejected!r}")
    h2=op.occupation(q2) if fin else None
    if fin:
        assert np.all((0.<=h2)&(h2<=1.))
        op.old=h0
    def invariants(U,h):
        M=U.reshape(2,op.n,2)
        counts=M[:,:,0].sum(axis=1)
        e=sum((c.mass+c.low)@M[j,:,0]+c.dk@M[j,:,1] for j,c in enumerate(op.k.c))
        if fin:
            counts+=op.measure@h
            e+=(op.measure*op.k.E)@h
        return np.r_[counts,e]
    inv0=invariants(U0,h0)
    errors=[float(max(abs(invariants(op.moment(q)[0],parent)/inv0-1.)))
            for q,parent in ((q1,h1),(q2,h2))]
    assert max(errors)<3e-10
    return q2,h2,dict(residual=max(r1["residual"],r2["residual"]),stage_1=r1,stage_2=r2,
                   max_relative_invariant_defect_at_each_stage=errors,
                   stage_2_seed=label,rejected_stage_2_seeds=rejected,
                   affine_daughter_target_outside_simple_cone_cells=nonreal,
                   minimum_and_maximum_affine_parent_target=affine_bounds,
                   first_stage_parent_min_max=None if h1 is None else [float(h1.min()),float(h1.max())])
