"""Number-weighted ballistic streaming of a thermal real scalar.
pip install numpy scipy; python thermal_ballistic.py
NOT a linear Boltzmann transfer-function calculation."""
import numpy as np
from scipy.constants import Boltzmann, e, c
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import zeta

def compute(m_eV=86.,gstar=55., H0=67.4,Om=.315,Or=9.2e-5, nq=1500,na=2400):
    Tnu=1.95*Boltzmann/e
    Tb=Tnu*(10.75/gstar)**(1/3)
    q=np.geomspace(1e-6,100,nq)
    w=q*q/np.expm1(q)
    N=trapezoid(w,q)
    assert abs(N/(2*zeta(3))-1)<1e-4
    a=np.geomspace(1e-12,1.,na)
    vb=[]
    for sl in np.array_split(np.arange(na),6):
        qq=q[None,:]
        aa=a[sl,None]*m_eV/Tb
        beta=qq/np.sqrt(qq*qq+aa*aa)
        vb.extend(trapezoid(beta*w[None,:],q,axis=1)/N)
    vb=np.array(vb)
    dist=np.r_[0.,cumulative_trapezoid((c/1000/H0)*vb/np.sqrt(Or+Om*a+(1-Om-Or)*a**4),a)]
    def report(z):
        az=1/(1+z)
        return (float(np.interp(az,a,dist)),float(np.interp(az,a,vb)*c/1000))
    return dict(m_eV=m_eV,gstar=gstar,Omega_h2=(2/3)*(m_eV/93.14)*(10.75/gstar),
                Delta_Neff=(4/7)*(10.75/gstar)**(4/3),
                mean_q_over_T=float(trapezoid(q*w,q)/N),
                z_eq=report(3400),z_recomb=report(1100),z_5=report(5),z_0=report(0))

if __name__=='__main__':
    for m in (29.,50.,86.,152.):
        x=compute(m,gstar=55*m/86)
        print(m,'eV',x)
    x=compute()
    y=compute(nq=900,na=1300)
    assert abs(x['z_eq'][0]/y['z_eq'][0]-1)<.005
    assert abs(x['mean_q_over_T']-2.70117803)<5e-5
    print('PASS: convergence and Bose momentum moment')