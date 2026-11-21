"""C1-smooth positive quadratically lowered Bose halo; η=ψ0=4 example.
Requires scipy numpy, sibling critical_lowered.py, lowered_virial.py.
DF f2=B(E)-B(E_t)-B'(E_t)(E-E_t) >=0 by strict convexity.
NOT an observed galaxy fit or proven globally stable system.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import PchipInterpolator
from critical_lowered import I as I1
from lowered_virial import J as J1

def I2(psi,eta=4.,nq=256):
    p=np.asarray(psi,dtype=float)
    B=1/np.expm1(eta)
    a=B*(1+B)
    val=I1(p,eta,nq)-2*a*p**2.5/15
    leading=4*a*(1+2*B)*p**3.5/105
    out=np.where(p<1e-4,leading,val)
    return float(out) if p.ndim==0 else out

def solve2(eta=4.,nq=256):
    ps=np.unique(np.r_[np.linspace(0,eta,450),
      eta-np.geomspace(1e-13,eta,420),eta])
    vals=I2(ps,eta,nq)
    assert abs(vals[0])<1e-20 and np.all(np.diff(vals)>0)
    Ic=float(vals[-1])
    eos=PchipInterpolator(ps,vals/Ic,extrapolate=False)
    eps=1e-5
    def rhs(r,w):return[-w[1]/r**2,r*r*float(eos(np.clip(w[0],0,eta)))]
    def edge(r,w):return w[0]
    edge.terminal=True;edge.direction=-1
    sol=solve_ivp(rhs,[eps,300],[eta-eps*eps/6,eps**3/3],
      events=edge,method='DOP853',rtol=1e-10,atol=1e-12,max_step=.04,dense_output=True)
    assert sol.success and len(sol.t_events[0])==1
    rt=float(sol.t_events[0][0])
    return Ic,rt,float(sol.sol(rt)[1]),sol,eos

def run():
    eta=4;Ic,rt,Q,sol,eos=solve2(eta)
    B=1/np.expm1(eta);a=B*(1+B)
    slope=np.pi/(2*np.sqrt(6)*Ic)
    xs=np.array([.002,.001,.0005])
    slopes=(1-eos(sol.sol(xs)[0]))/xs
    eps=np.geomspace(.001,.01,15)
    edge_exponent=float(np.polyfit(np.log(eps),np.log(I2(eps,eta)),1)[0])
    d=np.geomspace(1e-5,eta-.1,60)
    f2=1/np.expm1(eta-d)-B-a*d
    assert np.all(f2>0)
    assert abs(edge_exponent-3.5)<.02
    assert abs(slopes[1]/slope-1)<.002
    rr=np.r_[0,np.geomspace(1e-5,.5,400),np.linspace(.5,rt,6500)[1:]]
    psi=np.clip(sol.sol(np.maximum(rr,1e-5))[0],0,eta)
    q=np.maximum(sol.sol(np.maximum(rr,1e-5))[1],0)
    psi[0]=eta;q[0]=0
    rho=I2(psi,eta,nq=240)/Ic
    J2=J1(psi,eta,nq=240)-2*a/35*psi**3.5
    K=simpson(rr*rr*J2/Ic,x=rr)
    W=-simpson(rr*rho*q,x=rr)
    virial=(2*K+W)/(-W)
    assert abs(virial)<1e-5
    return dict(eta=eta,I2_c=Ic,nu2_c=4*Ic/np.sqrt(np.pi),
        central_slope=slope,measured_slopes=list(map(float,slopes)),
        r_t_over_R=rt,Q_total=Q,edge_exponent=edge_exponent,
        K=K,W=W,relative_virial_residual=virial,
        tests_passed=4,
        caveat="C1 at escape is not proof of Landau damping or stability.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
