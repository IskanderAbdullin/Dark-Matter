"""Virial and isotropic pressure identity for a critical lowered Bose halo.
Run in same directory as critical_lowered.py: python lowered_virial.py.
Requires numpy, scipy. No GitHub Actions and no empirical fitting.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import simpson
from critical_lowered import I,solve

def J(psi,eta=4.,nq=240):
    v=np.asarray(psi,dtype=float)
    if np.any(v<0) or np.any(v>eta):raise ValueError("0<=psi<=eta")
    nodes,weights=leggauss(nq)
    u=(nodes+1)/2; w=weights/2
    a=eta-v[...,None]*(1-u*u)
    d=v[...,None]*(1-u*u)
    diff=np.exp(a)*np.expm1(d)/(np.expm1(a)*np.expm1(eta))
    out=v**2.5*np.sum(w*u**4*diff,axis=-1)
    return float(out) if v.ndim==0 else out

def run():
    Ic,rt,Qt,res,eos=solve(4,nq=256)
    rr=np.r_[0.,np.geomspace(1e-5,.5,400),np.linspace(.5,rt,4800)[1:]]
    psi=np.clip(res.sol(np.maximum(rr,1e-5))[0],0.,4.)
    Q=np.maximum(res.sol(np.maximum(rr,1e-5))[1],0.)
    psi[0]=4;Q[0]=0
    x=I(psi,4,nq=240)/Ic
    Pbar=(2./3)*J(psi,4,nq=240)/Ic
    Kbar=simpson(rr**2*(1.5*Pbar),x=rr)
    Wbar=-simpson(rr*x*Q,x=rr)
    probe=np.array([.01,.1,1.,5.,15.,30.,rt*.99])
    ps=res.sol(probe)[0]
    hh=1e-5*np.minimum(1.,4-ps)
    hh=np.minimum(hh,ps*.4)
    def P(pp):return (2./3)*J(pp,4,nq=240)/Ic
    dPdpsi=(P(ps+hh)-P(ps-hh))/(2*hh)
    exact_x=I(ps,4,nq=240)/Ic
    err=float(np.max(np.abs(dPdpsi-exact_x)/exact_x))
    resid=float((2*Kbar+Wbar)/(-Wbar))
    assert abs(resid)<5e-5
    assert err<3e-6
    assert Pbar[-1]<1e-12
    return dict(eta=4,psi0=4,rt_over_R=rt,Q_total=Qt,
     dimensionless_kinetic_K=float(Kbar),
     dimensionless_potential_W=float(Wbar),
     virial_2K_plus_W_div_abs_W=resid,
     max_local_dP_dpsi_over_density_relative_error=err,
     pressure_edge_normalized=float(Pbar[-1]),
     checks_passed=3,
     warning="Isotropic DF virial/Jeans consistency, NOT stability proof or formation mechanism.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
