"""Support-edge-aligned physical Bose DG-P1 *moment representability* audit.

For a finite-parent cutoff P=2 GeV, the original first scalar cell crosses
the exact φ on-shell support floor. Split the cell exactly at the kinematic
threshold, keep the lower subcell unchanged and test the positive fixed-FD
scalar-only strong evolution in the active upper subcell. This is not a fully
coupled physical kinetic solver or a pointwise spectral accuracy proof.
"""
import math
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import root
from scipy.integrate import quad
from kinetics.phi_bose_scalar_fixed_bath_strong import mode_coefficients
from kinetics.phi_bose_moment_realisability_boundary import boundary_crossing
from kinetics.phi_bose_physical_cutcell_weak import MPHI,shell,cutcell_weak,scalar_edges


def tangent_segment(L,U,R,Q):
    """Physical two-hat susceptibility/tangent on a SINGLE arbitrary energy segment."""
    assert MPHI<=L<U
    C=np.array([(U*R-Q)/(U-L),(Q-L*R)/(U-L)])
    lo=math.sqrt(max(0,(L-MPHI)*(L+MPHI)))
    hi=math.sqrt((U-MPHI)*(U+MPHI))
    def integrand(p,i,j):
        E=math.hypot(p,MPHI)
        hR=(p*p/(E+MPHI)-(L-MPHI))/(U-L)
        hL=1-hR
        h=(hL,hR)
        b=1/math.expm1(E)
        return p*p/(2*math.pi**2)*h[i]*h[j]*b*(1+b)
    G=np.array([[quad(lambda p:integrand(p,i,j),lo,hi,
                      epsabs=1e-18,epsrel=3e-11)[0]
                 for j in (0,1)] for i in (0,1)])
    return C,G,np.linalg.solve(G,-C)


def positive_eta_moments(L,U,eta_left,eta_right,n=128):
    """Actual physical DOS two moments of positive affine inverse fugacity."""
    assert MPHI<=L<U and eta_left>0 and eta_right>0
    x,w=leggauss(n)
    k=L+(U-L)*(x+1)/2
    d=(U-L)*w/2
    hR=(k-L)/(U-L)
    eta=eta_left*(1-hR)+eta_right*hR
    B=1/np.expm1(eta)
    DOS=k*np.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
    return np.array([sum(d*DOS*B*(1-hR)),sum(d*DOS*B*hR)])


def target_moments(L,U,t,P=2.,nk=96,np_parent=48):
    """Actual positive scalar fixed-FD modal time moments, not frozen-source ray."""
    x,w=leggauss(nk)
    k=L+(U-L)*(x+1)/2
    d=(U-L)*w/2
    hR=(k-L)/(U-L)
    B0=1/np.expm1(k)
    delta=[]
    for ki in k:
        S,D,*_=mode_coefficients(float(ki),P,np_parent)
        z=-D*t
        exprel=math.expm1(z)/z if z!=0 else 1.
        delta.append(t*(S-D/math.expm1(ki))*exprel)
    db=np.array(delta)
    DOS=k*np.sqrt((k-MPHI)*(k+MPHI))/(2*math.pi**2)
    return np.array([sum(d*DOS*(B0+db)*(1-hR)),
                     sum(d*DOS*(B0+db)*hR)])


def fit_positive_etas(L,U,t,P=2.):
    """Construct exact two-moment positive affine-eta fit on the active subcell."""
    target=target_moments(L,U,t,P)
    baseline=positive_eta_moments(L,U,L,U)
    def fun(x):
        etaL=L*np.exp(x[0])
        etaR=U*np.exp(x[1])
        got=positive_eta_moments(L,U,etaL,etaR)
        return (got-target)/baseline
    sol=root(fun,[0.,0.],tol=1e-11)
    etaL=float(L*np.exp(sol.x[0]))
    etaR=float(U*np.exp(sol.x[1]))
    return dict(t=t,converged=bool(sol.success),
                fitted_eta_left=etaL,fitted_eta_right=etaR,
                relative_eta_shifts=[etaL/L-1,etaR/U-1],
                relative_moment_residual=fun(sol.x).tolist())


def audit():
    P=2.
    kfloor=shell(P)[1]
    out=[]
    for n in (24,32,40):
        val=cutcell_weak(P,10,16,n)
        m,U=val['edges'][:2]
        R,Q=val['bycell'][0,:2]
        C0,G0,old=tangent_segment(m,U,R,Q)
        C1,G1,new=tangent_segment(kfloor,U,R,Q)
        tstar=boundary_crossing(n)['boundary_time']
        out.append(dict(n=n,kinematic_floor=kfloor,original_edge=U,
                        Rphi=R,Qphi=Q,
                        original_eta_tangent=old.tolist(),
                        support_aligned_eta_tangent=new.tolist(),
                        original_frozen_ray_boundary_time=tstar,
                        positive_active_cell_fits=[
                            fit_positive_etas(kfloor,U,tstar*mult)
                            for mult in (1.1,2.,10.,100.)]))
    return dict(schema_version=1,
                scope="physical support-aligned first-bin positive fixed-fermion modal moment closure only, not full coupled strong",
                rows=out,class_ready=False)


if __name__=='__main__':
    print(json.dumps(audit(),indent=2,allow_nan=False))
