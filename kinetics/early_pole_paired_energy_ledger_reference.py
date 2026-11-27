"""Event-resolved parent energy demand and bin-centroid energy defect."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre, expit
from early_pole_paired_hartree_matrix_reference import DirectDriver
from early_pole_decay_mode_clock_reference import X,G,KAPPA


def ledger(d,x,u,b):
    m=d.gap(x,u)
    p=d.p
    E=np.hypot(p,m)
    A0=(x*x+m*m)/2
    D=(x*x-m*m)/2
    s=A0*E+p*D
    lo=(x*x*E*E+D*D)/s
    span=2*p*D/(m*m)
    ll=np.maximum(lo[:,None],E[:,None]+d.ef[None,:-1])
    hh=np.minimum((lo+span)[:,None],E[:,None]+d.ef[None,1:])
    dr=np.maximum(hh-ll,0.)
    assert x>m
    row=KAPPA*G*x*(x*x-m*m)/(E*p)
    A=np.log1p(np.exp(-ll)*(-np.expm1(-dr))/(1+np.exp(-ll-dr)))
    # Use quadrature only for the first moment inside each exact overlap.
    g,w=roots_legendre(6)
    J=np.zeros_like(ll)
    C=ll+.5*dr
    for z,wx in zip(g,w,strict=True):
        v=C+dr*z/2
        J+=wx/2*v*expit(-v)
    J*=dr
    mid=(d.fed[1:]+d.fed[:-1])/2
    # True energy-weighted events use the actual E_parent at each event.
    # Number closure has identical shared rates in both calculations.
    net=A*(1-b)+u[:,None]*(A-dr*b)
    q=d.w*row
    true=q@(J*(1-b)+u[:,None]*(J-dr*C*b)).sum(axis=1)
    eq=q@(net*(E[:,None]+mid)).sum(axis=1)
    gross=q@((J*(1-b)*(1+u[:,None]))).sum(axis=1)
    n=d.w*u
    nf=d.wf*b
    U=(m*m-(x*1e-9)**2)**2/(2*d.lam)
    rho=n@E-U+nf@mid
    trace=(x*1e-9)**2*(n@(1/E))
    return np.array([true,eq,gross,trace,rho])


def one(ns,nf):
    d=DirectDriver(ns,nf,7.65175e-8,0.,1e-9)
    T=np.log(10.)
    atol=np.r_[np.full(d.n,5e-8),np.full(len(d.wf),6e-12)]
    sol=solve_ivp(d.rhs,(0.,T),np.zeros(len(atol)),method='BDF',jac=d.jac,
                  atol=atol,rtol=1e-8,max_step=.020,dense_output=True)
    assert sol.success
    checks=[]
    for nr in (4,8):
        x,w=roots_legendre(nr)
        # The stiff solver's own adapted step boundaries resolve the ignition.
        h=(sol.t[1:]-sol.t[:-1])/2
        t=sol.t[:-1,None]+h[:,None]*(x[None,:]+1)
        ww=(h[:,None]*w).ravel()
        integ=[]
        for z in t.ravel():
            q=sol.sol(z); integ.append(ledger(d,X*np.exp(z),q[:d.n],q[d.n:]))
        integ=np.array(integ)
        V=ww@integ
        checks.append(V.tolist())
    true,pr,gross,trace,_=checks[-1]
    q=sol.y[:,-1]
    end=ledger(d,.05,q[:d.n],q[d.n:])
    assert np.max(abs(np.array(checks[-1])/checks[0]-1))<.003
    assert abs((pr+trace)/end[4]-1)<3e-5,checks
    return dict(scalar_bins=d.n,fermion_bins=len(d.wf),
                integrated_parent_net_energy_over_T4=true,
                integrated_centroid_energy_over_T4=pr,
                integrated_parent_gross_forward_energy_over_T4=gross,
                explicit_vacuum_trace_work_over_T4=trace,final_daughter_energy_over_T4=end[4],
                signed_projected_energy_relative_defect=pr/true-1,
                projected_continuity_defect=(pr+trace)/end[4]-1,
                multiquadrature_repeat=checks)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=65)
    p.add_argument('--fermion',type=int,default=90)
    q=p.parse_args()
    print(json.dumps(one(q.scalar,q.fermion),indent=2,allow_nan=False))
