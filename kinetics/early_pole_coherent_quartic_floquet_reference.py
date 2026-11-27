"""Conditional uniform real-field self-resonance from a stored scalar energy."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk


def edges(k2):
    return (3*(1-k2),1-2*k2+2*np.sqrt(1-k2+k2*k2))


def monodromy(k2,k_sq,rtol):
    # Background phi/A=cn(z|k2), z=Omega*t. Its squared period is
    # 2K, and the fluctuation responds to V'' in that period.
    P=2*ellipk(k2)
    n=len(k_sq)
    def rhs(z,y):
        a,b,c,d=y.reshape(4,n)
        cn=ellipj(z,k2)[1]
        L=k_sq+1-2*k2+6*k2*cn*cn
        return np.r_[c,d,-L*a,-L*b]
    y0=np.r_[np.ones(n),np.zeros(2*n),np.ones(n)]
    sol=solve_ivp(rhs,(0,P),y0,method='DOP853',rtol=rtol,atol=rtol*.1,max_step=.08)
    assert sol.success
    a,b,c,d=sol.y[:,-1].reshape(4,n)
    det=a*d-b*c
    tr=(a+d)*.5
    mu=np.arccosh(np.maximum(abs(tr),1.))/P
    return mu,float(np.max(abs(det-1))),P


def background(rho,m0,lam):
    # Stable positive root for rho=m0^2 A^2/2+lambda A^4/24.
    a=lam/24; b=m0*m0/2
    A2=2*rho/(b+np.sqrt(b*b+4*a*rho))
    O=np.sqrt(m0*m0+lam*A2/6)
    return np.sqrt(A2),O,lam*A2/(12*O*O)


def one(spec,M):
    r=spec['checkpoints'][-1]
    x=spec['x_end']
    rho=r['integrated_scalar_energy_over_T4']
    m0=1e-9*x
    A,O,k=background(rho,m0,7.65175e-8)
    l,h=edges(k)
    q2=np.linspace(l+.000002,h-.000002,121)
    compare=np.r_[q2,l-.003,h+.003,0.,.01,.4,3.]
    mu,det,P=monodromy(k,compare,2e-11)
    mu2,det2,_=monodromy(k,compare,1e-12)
    assert mu2[:len(q2)].max()>.03 and np.max(mu2[-6:])<3e-6
    assert abs(mu2.max()/mu.max()-1)<1e-7
    assert det2<2e-10
    ids=np.argmax(mu2)
    H=np.sqrt(np.pi**2*106.75/90)*(M/x)/2.435e18
    return dict(x=x,mass_ratio=spec['mass_ratio'],scalar_energy_over_T4=rho,
                assumed_homogeneous_amplitude_over_T=A,oscillator_Omega_over_T=O,
                elliptic_parameter=k,modulation_period_in_Omega_time=P,
                unstable_q_interval=[np.sqrt(l)*O,np.sqrt(h)*O],
                peak_q_over_T=float(np.sqrt(q2[ids])*O),
                maximum_floquet_exponent_over_Omega=float(mu2[ids]),
                conditional_growth_rate_over_radiation_reference_H=float(mu2[ids]*O/H),
                max_monodromy_determinant_relative_error=det2,
                max_growth_relative_change_across_tolerances=float(abs(mu2.max()/mu.max()-1)),
                check_outside_band_max_exponent=float(max(mu2[-6:])))


def check():
    A,O,k=background(.003,0.,.1)
    assert k==.5 or abs(k-.5)<1e-15
    a,b=edges(k)
    assert abs(a-1.5)<1e-14 and abs(b-np.sqrt(3))<1e-14
    mu,d,p=monodromy(k,np.linspace(a,b,160),1e-12)
    assert .035<mu.max()<.037
    return dict(known_pure_quartic_kappa_sq_band=[a,b],max_numerical_exponent=float(mu.max()),wronskian_error=d)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    p.add_argument('--parent-mass-GeV',type=float,default=86.)
    a=p.parse_args()
    out=dict(known_pure_quartic_check=check())
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        out[f.name]=one(json.loads(raw),a.parent_mass_GeV)
    print(json.dumps(out,indent=2,allow_nan=False))
