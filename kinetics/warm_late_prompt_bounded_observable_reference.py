"""Compare whole normalized phase distributions, beyond their shell moments."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_prompt_two_moment_evolution_reference import kin
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def phase(p,A,B,mass,g,order):
    # At fixed shell, the sign of f_A-f_B is the sign of an
    # affine entropy difference. Split before integrating its L1.
    k=kin(p,mass)
    rk=np.diff(k)
    a,b=A-B
    split=np.clip(np.divide(-a,b,out=np.zeros_like(a),where=b!=0.),0.,1.)
    lo=np.stack((np.zeros_like(a),split),axis=1)
    hi=np.stack((split,np.ones_like(a)),axis=1)
    x,w=roots_legendre(order)
    kval=lambda u:k[:-1,None]+rk[:,None]*u
    edge=lambda u:np.sqrt(kval(u)*(kval(u)+2*mass))
    ll,rr=edge(lo),edge(hi)
    step=(rr-ll)/2
    q=ll[:,:,None]+step[:,:,None]*(1+x)
    K=kin(q,mass)
    u=(K-k[:-1,None,None])/rk[:,None,None]
    TA=A[0,:,None,None]+A[1,:,None,None]*u
    TB=B[0,:,None,None]+B[1,:,None,None]*u
    assert np.min(TA)>0. and np.min(TB)>0.
    diff=a[:,None,None]+b[:,None,None]*u
    T0=np.minimum(TA,TB)
    qa,qb=np.exp(-TA),np.exp(-TB)
    da=np.exp(-T0)*(-np.expm1(-abs(diff)))
    sign=1 if g==1 else -1
    denom=(1-sign*qa)*(1-sign*qb)
    F=da/denom
    E=K+mass
    measure=g*q*q*step[:,:,None]*w/(2*np.pi**2)
    dx=measure*F
    v=q/E
    eps=np.array([.3,1.,3.,10.,30.,100.,300.])
    ang=np.sinc(v[...,None]*eps/np.pi)
    # Stable signed particle-measure change, evaluated together.
    signed=-dx*np.sign(diff)
    G=measure*np.exp(-TB)/(1-sign*np.exp(-TB))
    return float(dx.sum()),float(signed.sum()),(signed[...,None]*ang).sum(axis=(0,1,2)),(G[...,None]*ang).sum(axis=(0,1,2))


def upper(p,A,B,d,g,n):
    # Do not rely on quadrature for the positive L1 majorant. On
    # each short subpanel use the minimum possible theta and the
    # affine maximum |Delta theta| in a mean-value bound.
    k=kin(p,d)
    x=np.linspace(0.,1.,n+1)
    E=k[:-1,None]+np.diff(k)[:,None]*x
    mom=np.sqrt(E*(E+2*d))
    a=A[0,:,None]+A[1,:,None]*x
    b=B[0,:,None]+B[1,:,None]*x
    assert np.min(a)>0. and np.min(b)>0.
    delta=abs((A[0]-B[0])[:,None]+(A[1]-B[1])[:,None]*x)
    angle=np.maximum(delta[:,:-1],delta[:,1:])
    v=np.minimum(a,b)
    low=np.minimum(v[:,:-1],v[:,1:])
    B0=np.exp(-low)
    f=B0/(-np.expm1(-low)) if g==1 else B0/(1+B0)
    bound=np.minimum(f,angle*f*(1+(1 if g==1 else -1)*f))
    l,r=mom[:,:-1],mom[:,1:]
    volume=g/(6*np.pi**2)*(r-l)*(r*r+r*l+l*l)
    return float(np.sum(volume*bound))


def run(active,vacuum):
    A,B=map(read,(active,vacuum))
    X,Y=A["snapshot"],B["snapshot"]
    assert A["hash"]==digest(X) and B["hash"]==digest(Y)
    assert A["time"]==B["time"] and B["plan"]["matching_active_grid_sha256"]==A["plan"]["seed_sha256"]
    p=np.array(X["lab_momentum_edges_GeV"])
    np.testing.assert_array_equal(p,Y["lab_momentum_edges_GeV"])
    keys=("entropy_variable_scalar","entropy_variable_fermion")
    a,b=([np.asarray(S[k]) for k in keys] for S in (X,Y))
    N0,N1=(get_mom(p,z)[0].sum(axis=1) for z in (a,b))
    out=[]
    for j,(d,g) in enumerate(((mp,1),(m,4))):
        states=[phase(p,a[j],b[j],d,g,o) for o in (12,24,48)]
        assert abs(states[-1][0]/states[-2][0]-1.)<1e-6
        q,net,dchar,char=states[-1]
        N,old=N0[j],N1[j]
        assert abs(net-(N-old))/old<4e-14
        delta=dchar/N-char*net/(N*old)
        cap=[upper(p,a[j],b[j],d,g,n) for n in (8,32,128)]
        assert q<=min(cap)*(1+1e-8)
        norm=2*cap[-1]/N
        assert np.max(abs(delta))<=norm*(1+1e-13)
        out.append(dict(estimated_L1_before_number_renormalization=q/old,
                        certified_formula_upper_before_number_renormalization=[x/old for x in cap],
                        normalized_L1_numerical_majorant=norm,
                        sampled_isotropic_velocity_characteristic_changes=delta.tolist()))
    return dict(scope="Phase-measure comparison on the same reconstructed finite mesh. The numerical positive "
                      "majorant bounds bounded single-particle test functions up to floating evaluation; "
                      "the sampled instantaneous velocity factor is not a cosmological transfer function.",
                source_sha256=[A["hash"],B["hash"]],time_H=A["time"],hard_GeV=A["plan"]["hard_GeV"],
                velocity_frequencies=[.3,1.,3.,10.,30.,100.,300.],species=["scalar","light"],controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("active",type=Path);p.add_argument("vacuum",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.active,a.vacuum),indent=2,allow_nan=False))
