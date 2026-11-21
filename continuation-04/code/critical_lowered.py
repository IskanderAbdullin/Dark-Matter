"""Exactly critical lowered-Bose finite-mass equilibrium and convex energy-Casimir checks.
Not a nonlinear stability theorem, thermodynamic equilibrium, or observational fit.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

def I(psi,eta=4.,nq=512):
    v=np.asarray(psi,dtype=float)
    if np.any((v<0)|(v>eta)):raise ValueError("0<=psi<=eta")
    nodes,ws=leggauss(nq)
    u=(nodes+1)/2
    w=ws/2
    a=eta-v[...,None]*(1-u*u)
    d=v[...,None]*(1-u*u)
    difference=np.exp(a)*np.expm1(d)/(np.expm1(a)*np.expm1(eta))
    integral=v**1.5*np.sum(w*u*u*difference,axis=-1)
    return float(integral) if v.ndim==0 else integral

def casimir(f,eta=4.):
    f=np.asarray(f,dtype=float)
    c=1/np.expm1(eta)
    def H(x):return (x+1)*np.log(x+1)-x*np.log(x)
    return eta*f-(H(f+c)-H(c))

def casimir_second(f,eta=4.):
    c=1/np.expm1(eta)
    return 1/((f+c)*(1+f+c))

def solve(eta=4.,minimum_y=1e-14,nq=512):
    psi=np.unique(np.r_[np.linspace(0,eta,450),
      eta-np.geomspace(minimum_y,eta,420),eta])
    vals=I(psi,eta,nq)
    assert vals[0]==0 and np.all(np.diff(vals)>0)
    Ic=float(vals[-1])
    eos=PchipInterpolator(psi,vals/Ic,extrapolate=False)
    eps=1e-5
    def rhs(xi,state):
        psi,Q=state
        return[-Q/xi**2,xi**2*float(eos(np.clip(psi,0,eta)))]
    def edge(xi,state):return state[0]
    edge.terminal=True
    edge.direction=-1
    res=solve_ivp(rhs,[eps,200],[eta-eps**2/6,eps**3/3],
       events=edge,dense_output=True,rtol=8e-10,atol=2e-11,max_step=.08)
    if not res.success or not len(res.t_events[0]):
        raise RuntimeError("Could not locate finite escape boundary")
    rt=float(res.t_events[0][0])
    return Ic,rt,float(res.sol(rt)[1]),res,eos

def run():
    eta=4.
    Ic,rt,Qt,res,eos=solve(eta)
    slope0=np.pi/(2*np.sqrt(6)*Ic)
    xs=np.array([.002,.001,.0005])
    y,q=res.sol(xs)
    slopes=(1-eos(y))/xs
    small_y=np.array([1e-6,1e-7,1e-8])
    coefficients=(Ic-I(eta-small_y,eta))/np.sqrt(small_y)
    assert abs(coefficients[-1]-np.pi/2)<1e-3
    assert abs(slopes[-1]/slope0-1)<.002
    assert abs(Ic-1.087257505)<1e-6
    c=1/np.expm1(eta)
    fmax=1/np.expm1(.2)-c
    f=1.32;h=1e-5
    cp=lambda f: eta-np.log1p(1/(f+c))
    derivative=(cp(f+h)-cp(f-h))/(2*h)
    assert abs(derivative/casimir_second(f,eta)-1)<1e-9
    assert casimir(fmax,eta)>0
    return dict(eta=eta,psi0=eta,I_critical=Ic,
      critical_density_slope=slope0,
      radial_slopes=[dict(xi=float(x),estimate=float(s)) for x,s in zip(xs,slopes)],
      sqrt_y_coeff=[dict(y=float(v),estimate=float(c)) for v,c in zip(small_y,coefficients)],
      boundary_rt_over_R=rt,total_dimensionless_mass_Q=Qt,
      casimir=dict(c=c,fmax_subcritical_eta4_psi03_8=fmax,
        C_double_prime_at_f0=float(casimir_second(0,eta)),
        C_double_prime_at_fmax=float(casimir_second(fmax,eta)),
        C_double_prime_finite_diff=float(derivative)),
      tests_passed=5,
      warning="Exactly critical ideal zero-momentum singularity integrable in density; no nonlinear Vlasov or collisional stability claim.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
