"""Critical Hartree-Fock Bose halo toy solver; run: pip install numpy scipy mpmath."""
import numpy as np
from mpmath import mp
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.constants import G, hbar, c, e, parsec

mp.dps = 35
NU=float(mp.zeta(1.5))
C=2*np.sqrt(np.pi)/NU
def equation_of_state(b):
    t=np.r_[0.,np.geomspace(1e-17,50.,850)]
    x=np.array([float(mp.polylog(1.5,mp.exp(-mp.mpf(str(z)))))/NU for z in t])
    y=t+b*(1-x)
    selected=np.r_[True,np.diff(y)>0]
    return PchipInterpolator(y[selected],x[selected],extrapolate=False)

def profile(b, rmax=8):
    eos=equation_of_state(b)
    eps=1e-6
    def rhs(r,w):
        y,q=w
        return [q/r**2, r*r*float(eos(max(0,y)))]
    s=solve_ivp(rhs,[eps,rmax],[eps**2/6,eps**3/3],rtol=1e-9,atol=1e-12,dense_output=True,max_step=.025)
    assert s.success,s.message
    return s,eos

def physical_benchmark():
    m=86*e/c**2; R=600*parsec; v=44.1e3; sm=0.125
    rho=v*v/(4*np.pi*G*R*R)
    a=np.sqrt(sm*m/(8*np.pi))
    g=4*np.pi*hbar*hbar*a/m
    b=2*g*rho/(m*m*v*v)
    return b,np.sqrt(6)*C*b, rho*1e-3

if __name__=='__main__':
    for b in (0,.001,.01,.05):
        solution,eos=profile(b)
        y,q=solution.sol(.03)
        print('b=',b,'rho(.03R)/rho0=',float(eos(y)),'u(.03R)=',float(np.sqrt(q/.03)))
    print('critical nu',NU,'ideal slope',C/np.sqrt(6))
    print('b, crossover xi, density g/cm3',physical_benchmark())
    assert abs(C/np.sqrt(6)-.553979545797)<1e-10
    s0,e0=profile(0); s1,e1=profile(.001)
    for r in (.5,1,2,4):
        assert abs(float(e0(s0.sol(r)[0]))-float(e1(s1.sol(r)[0])))<.003
    print('PASS: reference checks')