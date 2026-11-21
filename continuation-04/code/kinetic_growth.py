"""Growing-mode imaginary-axis Vlasov–Poisson dispersion for homogeneous BE gas.
No collisions, one noncondensed species. NOT a global finite halo stability test.
"""
import json,math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from mpmath import mp
mp.dps=40
nu=2.61
critical=float(mp.zeta(mp.mpf("1.5")))
lo,hi=mp.mpf(0),mp.mpf(10)
for _ in range(140):
    md=(lo+hi)/2
    if float(mp.polylog(mp.mpf("1.5"),mp.exp(-md)))>nu:lo=md
    else:hi=md
delta=float((lo+hi)/2)
I_norm=math.sqrt(math.pi)*nu/4

def occupation(q):return 1/math.expm1(delta+q*q)

def angular(q,s):
    x=q/s
    if x<.01:
        xx=x*x
        return xx/3-xx*xx/5+xx**3/7-xx**4/9
    return 1-math.atan(x)/x

def A_imag(s):
    if s==0:return float(mp.polylog(mp.mpf(".5"),mp.exp(-delta))/nu)
    def integrand(q):
        f=occupation(q)
        return q*q*f*(1+f)*angular(q,s)
    edges=[0,.0001,.001,.01,.1,.5,1,2,4,8,15]
    return sum(quad(integrand,edges[i],edges[i+1],
        epsabs=1e-8,epsrel=5e-7,limit=300)[0] for i in range(len(edges)-1))/I_norm

def growth(kR,b=0):
    # gamma_hat = gamma R/V0, s = gamma_hat/(sqrt(2)*kR)
    k2=kR*kR
    A0=A_imag(0)
    if k2>=A0/(1+b*A0):return 0.
    def eq(s):
        A=A_imag(s)
        return A/(1+b*A)-k2
    upper=1.
    while eq(upper)>0:upper*=2
    s=brentq(eq,1e-13,upper,xtol=1e-11,rtol=1e-11)
    return s*kR*math.sqrt(2)

def run():
    vals=[]
    for k in [.01,.1,1,3,5,10,15,20,25,30,31,31.7,32]:
        gam=growth(k)
        vals.append(dict(kR=k,wavelength_over_R=2*math.pi/k,
            gamma_R_over_V0=gam,
            e_fold_Myr_at_R06kpc_V044kms=13.302553/gam if gam>0 else None))
    assert abs(growth(.01)-1)<.0002
    assert abs(A_imag(1e-8)/A_imag(0)-1)<.005
    assert vals[-1]["gamma_R_over_V0"]==0
    return dict(nu=nu,delta=delta,static_A=A_imag(0),
      kR_threshold=math.sqrt(A_imag(0)),imaginary_growth=vals,
      tests_passed=3,
      limitations="Homogeneous infinite Vlasov medium with Jeans background subtraction. No global halo eigenmodes, baryons, Bose stars, galaxy fit or self-consistent perturbations.")

if __name__=="__main__":print(json.dumps(run(),indent=2))
