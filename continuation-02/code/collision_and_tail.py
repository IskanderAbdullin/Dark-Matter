"""Fiducial collision proxies + far-field critical Bose halo test.
pip install numpy scipy mpmath; python collision_and_tail.py
Depends on mean_field_regularization.py in this directory.
"""
import numpy as np
from scipy.integrate import trapezoid,solve_ivp
from scipy.constants import G,parsec
from mean_field_regularization import equation_of_state,profile

R=600*parsec
V0=44.1e3
sigma_over_m=.125 # m^2/kg = 1.25 cm^2/g
rho0=V0**2/(4*np.pi*G*R**2)
lambda0=1/(rho0*sigma_over_m)
t0=lambda0/V0/(365.25*86400*1e9)
s,eos=profile(0,rmax=60)
print('rho0 g/cm3',rho0/1000,'central mfp kpc',lambda0/(1000*parsec),
      'central collision time Gyr',t0)
for rr in (.1,1,3,6,10,20):
    x=float(eos(s.sol(rr)[0]))
    z=np.geomspace(1e-5,rr,2000)
    tau=rho0*sigma_over_m*R*trapezoid(eos(s.sol(z)[0]),z)
    print('xi, rho/rho0, coll time Gyr, radial tau',rr,x,t0/x,float(tau))
# solve in log-radius t=log(xi), p=xi*dy/dxi
eos=equation_of_state(0.)
def ode(t,w):
    y,p=w
    return [p,np.exp(2*t)*float(eos(y))-p]
eps=1e-5
res=solve_ivp(ode,[np.log(eps),np.log(1e5)],[eps**2/6,eps**2/3],
              dense_output=True,rtol=1e-8,atol=1e-10,max_step=.02)
assert res.success
for rr in (100,1000,10000,100000):
    y,p=res.sol(np.log(rr))
    print('tail xi, xi²rho/rho0, u²',rr,float(rr**2*eos(y)),float(p))
