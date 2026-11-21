"""Homogeneous local Vlasov–Poisson zero-frequency Jeans threshold for an ideal BE gas.
Requires Python 3, scipy, numpy, mpmath. NOT a global halo stability proof.
"""
import json
import mpmath as mp
import numpy as np
from scipy.constants import G,hbar,c,elementary_charge,parsec
from scipy.integrate import quad
mp.mp.dps=45
nuc=mp.zeta(mp.mpf("1.5"))
R=600*parsec;V=44.1e3
m=86*elementary_charge/c**2
rho=V**2/(4*np.pi*G*R**2)
sigma_per_m=.125
scattering_a=np.sqrt(sigma_per_m*m/(8*np.pi))
g=4*np.pi*hbar**2*scattering_a/m
b=2*g*rho/(m*m*V*V)

def t_from_nu(nu):
    target=mp.mpf(str(nu))
    if float(nu)>=float(nuc):return mp.mpf(0)
    lo,hi=mp.mpf(0),mp.mpf(30)
    for _ in range(180):
        md=(lo+hi)/2
        if mp.polylog(mp.mpf("1.5"),mp.exp(-md))>target:lo=md
        else:hi=md
    return (lo+hi)/2

def analytic_A(delta):
    z=mp.exp(-delta)
    return mp.polylog(mp.mpf(".5"),z)/mp.polylog(mp.mpf("1.5"),z)

def momentum_quad_A(delta):
    d=float(delta)
    def occ(q):return 1/np.expm1(d+q*q)
    breaks=sorted(set([0,max(np.sqrt(d)/10,1e-8),np.sqrt(d),
      10*np.sqrt(d),.01,.1,.5,2.,6.,12.]) | {0.0,12.0})
    breaks=[x for x in breaks if x>=0 and x<=12]
    n=density=0.
    for low,high in zip(breaks[:-1],breaks[1:]):
        density+=quad(lambda q:q*q*occ(q),low,high,epsabs=1e-10,epsrel=2e-10,limit=250)[0]
        n+=quad(lambda q:q*q*occ(q)*(1+occ(q)),low,high,
                epsabs=1e-9,epsrel=2e-10,limit=250)[0]
    return n/density

def run():
    rows=[]
    for nu in [.03,.5,1.,1.5,2.2,2.61,float(nuc)]:
        delta=t_from_nu(nu)
        A=analytic_A(delta) if delta>0 else mp.inf
        if A==mp.inf:
            length0=0.;lengthHF=float(2*np.pi*R*np.sqrt(b)/parsec)
        else:
            length0=float(2*mp.pi*R/mp.sqrt(A)/parsec)
            lengthHF=float(2*mp.pi*R*mp.sqrt(1/A+b)/parsec)
        rows.append(dict(nu=float(nu),delta=float(delta),z0=float(mp.exp(-delta)),
           A=None if A==mp.inf else float(A),
           kJ_times_R_ideal=None if A==mp.inf else float(mp.sqrt(A)),
           lambdaJ_ideal_pc=length0,lambdaJ_HF_pc=lengthHF))
    delta=t_from_nu(2.61)
    direct=momentum_quad_A(delta);analytic=float(analytic_A(delta))
    assert abs(direct/analytic-1)<1e-7
    assert abs(rows[-2]["lambdaJ_ideal_pc"]-118.46908658)<1e-5
    assert abs(rows[-1]["lambdaJ_HF_pc"]-.00182989369028)<1e-9
    assert abs(rows[0]["lambdaJ_ideal_pc"]-2*np.pi*600)<50
    return dict(source="Abdullin & Popov JCAP 2021 ideal Bose DF",
        background="infinite homogeneous Vlasov–Poisson Jeans subtraction",
        nu_c=float(nuc),b_repulsive_HF=b,R_kpc=.6,V0_km_s=44.1,
        direct_quad_A_nu_2_61=direct,polylog_A_nu_2_61=analytic,
        relative_quad_error=abs(direct/analytic-1),checks_passed=4,
        profiles=rows,
        caveat="Local homogenous marginal threshold, not the full halo's stability or a galaxy fit.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
