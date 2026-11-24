"""Physical upper-support Bose DG-P1 inverse-eta endpoint curvature gate.

The true pointwise eta_dot is positive (thermal net inverse absorption),
but physical DOS weighted AFFINE-eta projected tangent can be negative
at the last active daughter energy support edge, even if that cell is
aligned exactly with kmax. Negative scales as width² via physical
curvature, not evidence for creation, BEC or finite-time strong evolution.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from kinetics.phi_bose_scalar_fixed_bath_strong import mode_coefficients
from kinetics.phi_bose_physical_cutcell_weak import MPHI,shell,cutcell_weak
from kinetics.phi_bose_support_aligned_moment_closure import tangent_segment


def modal_eta_tangent(k,P=2.,nparent=48):
    """True scalar inverse fugacity rate of exact fixed FD1/FD2 mode."""
    S,D,*_=mode_coefficients(k,P,nparent)
    b=1/math.expm1(k)
    return -(S-D*b)/(b*(1+b))


def physical_weak_segment(L,U,P=2.,nk=48,nparent=48):
    """Independent scalar-first conditional weak number/energy source on segment."""
    z,w=leggauss(nk)
    ks=L+(U-L)*(z+1)/2
    dk=(U-L)*w/2
    R=Q=0.
    for k,wi in zip(ks,dk):
        S,D,*_=mode_coefficients(float(k),P,nparent)
        b=1/math.expm1(float(k))
        dos=k*math.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
        val=wi*dos*(S-D*b)
        R+=val;Q+=k*val
    return R,Q


def curvature_fit(P=2.):
    K=shell(P)[2]
    ds=np.array([.001,.002,.003,.004,.005,.006])
    # g(K−δ)=Aδ+Bδ²+..., so g/δ=A+Bδ+Cδ².
    y=np.array([modal_eta_tangent(K-d,P)/d for d in ds])
    c=np.polyfit(ds,y,2)
    return float(c[2]),float(c[1])


def audit(N=32,P=2.):
    v=cutcell_weak(P,10,16,N)
    K=shell(P)[2]
    j=next(j for j,(lo,hi) in enumerate(zip(v['edges'][:-1],v['edges'][1:]))
           if lo<K<=hi)
    base=v['edges'][j]
    A,B=curvature_fit(P)
    rows=[]
    for nsub in (1,2,4,8,16,32,64):
        h=(K-base)/nsub
        R,Q=physical_weak_segment(K-h,K,P,48,48)
        C,G,eta=tangent_segment(K-h,K,R,Q)
        predicted=-B*h*h/6
        rows.append(dict(subdivisions=nsub,subcell_width_GeV=h,
                         Rphi_GeV4=R,Qphi_GeV5=Q,
                         eta_left_dot_GeV=float(eta[0]),
                         eta_right_dot_GeV=float(eta[1]),
                         leading_curvature_prediction_right=predicted,
                         ratio_to_asymptotic=float(eta[1]/predicted)))
    return dict(schema_version=1,N=N,P=P,
                last_active_bin=j,kmax=K,original_edge=base,
                physical_modal_eta_dot_curvature_A=A,B=B,
                rows=rows,
                scope='physical upper support DG-P1 instantaneous negative only, NOT physical Bose growth/strong finite-time',
                class_ready=False)


if __name__=='__main__':
    print(json.dumps(audit(),indent=2,allow_nan=False))
