"""Critical six-dimensional Bose phase-volume exponents under finite Hartree-Fock b.
Requires numpy scipy mpmath. Run: python ir_phase_volume_exponents.py.
These are local asymptotic checks, not global galactic stability.
"""
import json
import numpy as np
import mpmath as mp
from scipy.special import beta,zeta
nu=float(zeta(1.5))
c=2*np.sqrt(np.pi)/nu
C0=4*np.sqrt(3)*np.pi**3
Cb0=32*np.sqrt(3)*np.pi**2*beta(.75,2.5)*c**1.5

def solve_t(y,b):
    mp.mp.dps=50
    yy=mp.mpf(str(y));bb=mp.mpf(str(b));nc=mp.zeta(mp.mpf("1.5"))
    def eq(s):return s*s+bb*(1-mp.polylog(mp.mpf("1.5"),mp.exp(-s*s))/nc)-yy
    lo=mp.mpf(0);hi=mp.sqrt(yy)
    for _ in range(200):
        md=(lo+hi)/2
        if eq(md)>0:hi=md
        else:lo=md
    return float(((lo+hi)/2)**2)

def run(b=.01):
    expected=1/(b*c)**2
    checks=[]
    for y in (1e-4,1e-5,1e-6,1e-7,1e-8):
        t=solve_t(y,b)
        checks.append(dict(y=y,t=t,ratio=(t/y**2)/expected))
    assert abs(checks[-1]["ratio"]-1)<.001
    L=np.log1p(1/np.array([1e3,1e4,1e5,1e6,1e7]))
    F=1/np.expm1(L)
    f1=C0*L**3;f2=Cb0*b**1.5*L**2.25
    s1=float(np.polyfit(np.log(F),np.log(f1),1)[0])
    s2=float(np.polyfit(np.log(F),np.log(f2),1)[0])
    assert abs(s1+3)<.01 and abs(s2+2.25)<.01
    return dict(b=b,critical_nu=nu,c=c,
      leading_ideal_volume_coefficient=C0,
      leading_HF_coefficient_for_b1=Cb0,
      implicit_EOS=checks,ideal_slope=s1,HF_slope=s2,
      Lp_ideal_threshold=3,Lp_HF_threshold=2.25,
      checks_passed=3,
      caveat="Fixed-b local Hartree-Fock limit F->infinity; very likely beyond microscopic validity at physical 86-eV coupling.")
if __name__=="__main__":print(json.dumps(run(),indent=2))
