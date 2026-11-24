"""Exact physical parent-split cut-cell Bose DG-P1 first-cell tangent control.

Uses independently integrated positive p²dp Bose susceptibility and physical
cut-cell weak Rphi/Qphi. Negative inverse-fugacity tangents are numerical
representability obstructions, NOT a condensate claim or finite-time solver.
"""
import math
import json
import numpy as np
from scipy.integrate import quad
from kinetics.phi_bose_physical_cutcell_weak import MPHI,cutcell_weak,scalar_edges


def first_cell_susceptibility(upper,temperature=1.,abs_tol=2e-18):
    """Positive 2x2 Gram matrix of local scalar inverse-eta moments."""
    m=MPHI
    pmax=math.sqrt(upper*upper-m*m)
    def g(p,i,j):
        e=math.hypot(p,m)
        hats=((upper-e)/(upper-m),(e-m)/(upper-m))
        bose=1/math.expm1(e/temperature)
        return p*p/(2*math.pi**2)*hats[i]*hats[j]*bose*(1+bose)
    cuts=[0.]
    for p in (m,10*m,100*m,1e-5,1e-4):
        if 0<p<pmax:cuts.append(p)
    cuts.append(pmax)
    cuts=sorted(set(cuts))
    G=np.zeros((2,2))
    for i in range(2):
        for j in range(i,2):
            G[i,j]=sum(quad(lambda p:g(p,i,j),a,b,
                            epsabs=abs_tol,epsrel=2e-11,limit=200)[0]
                       for a,b in zip(cuts[:-1],cuts[1:]))
            G[j,i]=G[i,j]
    return G


def tangent(n,pmax=2.,nrad=10,nangle=16):
    """Actual thermal first-bin source and dual moment tangent dη=-G^-1 C."""
    v=cutcell_weak(pmax,nrad,nangle,n)
    lo,up=v['edges'][:2]
    R,Q=v['bycell'][0,:2]
    C=np.array([(up*R-Q)/(up-lo),(Q-lo*R)/(up-lo)])
    G=first_cell_susceptibility(up)
    deta=np.linalg.solve(G,-C)
    bound=lo/(-deta[0]) if deta[0]<0 else None
    return dict(N=n,pmax=pmax,first_bin_top_GeV=float(up),
                Rphi_first_GeV4=float(R),Qphi_first_GeV5=float(Q),
                source_hat_left_right=C.tolist(),G=G.tolist(),
                eigenvalues_G=np.linalg.eigvalsh(G).tolist(),
                deta_left_right=deta.tolist(),
                linear_tangent_positive_eta_limit_GeV_inverse=bound,
                first_cell_kinematically_accessible=bool(R!=0.))


def audit():
    return dict(schema_version=1,scope="one-bin instantaneous projected tangent only, not physical nonlinear time or Bose condensation",
                rows=[tangent(n) for n in (8,12,16,20,24,32,40,48,64)],
                class_ready=False)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,allow_nan=False))
