"""A physical Bose Bernstein-P2 projected tangent can be pointwise negative.

For each physical first ACTIVE scalar cell [kmin(P), first quadratic edge],
integrate the true on-shell scalar fixed-fermion gain/loss in the independent
energy-first + inverse-parent order. The physical susceptibility G_ab uses
actual k*sqrt(k²-m²)dk/(2π²) and the nonnegative degree-two Bernstein basis.

Even when each exact Bose mode is DEPLETING and true eta_dot(k)>0, its
quadratic moment projection need not be pointwise positive. This is a
numerical spectral representation counterexample, NOT physical Bose creation
or condensate, nor a demonstration of a full-time solver failure.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss

from kinetics.phi_bose_physical_cutcell_weak import MPHI,shell,scalar_edges
from kinetics.phi_bose_invariant_parent_support import global_support_floor
from kinetics.phi_bose_scalar_fixed_bath_strong import mode_coefficients
from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled,bernstein


def first_active_cell(P,nphi):
    """Physical shell lower support, original first scalar quadratic bin top."""
    lower=global_support_floor(MPHI,P)
    upper=float(scalar_edges(P,nphi)[1])
    if upper<=lower:
        raise ValueError("Original first bin is physically inactive")
    return lower,upper


def independent_bernstein_p2_tangent(P,nphi,nk=64,nparent=64):
    """Independent k-first weak R/Q mapped into physical 3-moment eta tangent."""
    L,U=first_active_cell(P,nphi)
    z,w=leggauss(nk)
    k=L+(U-L)*(z+1)/2
    dk=(U-L)*w/2
    x=(k-L)/(U-L)
    B=bernstein(x)
    bose=1/np.expm1(k)
    dos=k*np.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
    W=dk*dos
    src=[]
    for kval,bval in zip(k,bose):
        S,D,*_=mode_coefficients(float(kval),P,nparent)
        src.append(S-D*bval)
    source=np.asarray(src)
    G=np.einsum("i,ia,ib->ab",W*bose*(1+bose),B,B)
    C=np.einsum("i,ia->a",W*source,B)
    deta=np.linalg.solve(G,-C)
    return dict(P=P,nphi=nphi,energy_interval=[L,U],
                d_eta_coefficients=deta.tolist(),
                physical_firstcell_Rphi=float(np.sum(W*source)),
                physical_firstcell_Qphi=float(np.dot(W*source,k)),
                min_eigenvalue_G=float(np.min(np.linalg.eigvalsh(G))),
                projected_minimum=projected_minimum(deta,L,U))


def projected_minimum(three_coefficients,L,U):
    """Exact interval minimum of a real quadratic Bernstein polynomial."""
    a,b,c=np.asarray(three_coefficients,float)
    curv=a-2*b+c
    choices=[(0.,a),(1.,c)]
    if curv>0:
        x=float(np.clip((a-b)/curv,0,1))
        q=a*(1-x)**2+2*b*x*(1-x)+c*x*x
        choices.append((x,q))
    x,val=min(choices,key=lambda p:p[1])
    return dict(value=float(val),x=float(x),E_GeV=float(L+x*(U-L)),
                sign="negative" if val<0 else "nonnegative")


def true_eta_dot(k,P,nparent=64):
    """Exact positive physical modewise inverse-fugacity rate at k."""
    S,D,*_=mode_coefficients(k,P,nparent)
    bose=1/math.expm1(k)
    return -(S-D*bose)/(bose*(1+bose))


def event_order_bernstein_tangent(P,nphi,nrad=8,nangle=16,nq=32):
    """Cross-check using independent common p-first three-species event operator."""
    model=Coupled(P=P,Nphi=nphi,nrad=nrad,nangle=nangle,nq=nq)
    source=model.source(model.initial_etas())["C"][2].reshape((-1,3))
    phi=model.species[2]
    G=phi.evaluate(model.initial_etas()[2],with_G=True)["G"]
    lo,hi=first_active_cell(P,nphi)
    j=next(i for i,(a,b) in enumerate(zip(phi.edges[:-1],phi.edges[1:]))
           if abs(a-lo)<1e-11 and abs(b-hi)<1e-11)
    deriv=np.linalg.solve(G[j],-source[j])
    return dict(d_eta_coefficients=deriv.tolist(),
                projected_minimum=projected_minimum(deriv,lo,hi),
                firstcell_C=source[j].tolist(),events=model.event_count)


def audit():
    rows=[]
    for P,nphi in ((2.,8),(2.,12),(8.,8),(8.,12),
                   (8.,16),(16.,16),(16.,32)):
        independently=independent_bernstein_p2_tangent(P,nphi,80,80)
        independently["event_order_8x16"]=event_order_bernstein_tangent(
            P,nphi,8,16,32)
        m=independently["projected_minimum"]
        independently["true_physical_eta_dot_at_projected_min"]=true_eta_dot(
            m["E_GeV"],P,80)
        rows.append(independently)
    return dict(schema_version=1,
                scope="single first physically active scalar cell, fixed FD thermal trial, instantaneous projection ONLY",
                rows=rows,class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
