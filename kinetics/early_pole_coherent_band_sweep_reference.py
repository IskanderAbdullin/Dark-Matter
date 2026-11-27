"""Conditional partial-coherence and narrow-band crossing diagnostic."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj, ellipk, roots_legendre

from early_pole_coherent_quartic_floquet_reference import background


def gap(k):
    # Exact narrow positive gap, rationalized to avoid cancellation.
    root=np.sqrt(1-k+k*k)
    lower=3*(1-k)
    width=3*k*k/(2*root+2-k)
    return lower,width


def profile(k):
    low,W=gap(k)
    z,w=roots_legendre(28)
    eta=(1+z)/2
    P=2*ellipk(k)
    if k<8e-4:
        mu=W*np.sqrt(eta*(1-eta))/4
        return float(W/8),float(.5*w@mu*W),'narrow_gap_asymptote'
    kap=low+W*eta
    q=len(kap)
    def rhs(t,y):
        a,b,c,d=y.reshape(4,q)
        L=kap+1-2*k+6*k*ellipj(t,k)[1]**2
        return np.r_[c,d,-L*a,-L*b]
    sol=solve_ivp(rhs,(0,P),np.r_[np.ones(q),np.zeros(q*2),np.ones(q)],
                  method='DOP853',rtol=6e-13,atol=4e-14,max_step=.04)
    a,b,c,d=sol.y[:,-1].reshape(4,q)
    det=a*d-b*c
    rad=((a-d)/2)**2+b*c
    assert min(rad)>-5e-14 and np.max(abs(det-1))<2e-11
    mu=np.arcsinh(np.sqrt(np.maximum(rad/det,0.)))/P
    return float(mu.max()),float(.5*w@mu*W),'numerical_monodromy'


def one(s):
    # An illustrative interpolation: the incoherent pole contribution
    # scales linearly with the fraction not assigned to the zero mode.
    r=s['checkpoints'][-1]; x=r['x']; lam=7.65175e-8
    rho=r['integrated_scalar_energy_over_T4']; pole=r['scalar_pole_over_T']
    m0=1e-9*x
    H=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
    g=r['net_log_growth_rate']
    rows=[]
    for f in (1e-6,1e-5,3e-5,1e-4,1e-3,.003,.01,.03,.10,.30,1.):
        bg=np.sqrt(m0*m0+(1-f)*(pole*pole-m0*m0))
        a,O,k=background(rho*f,bg,lam)
        low,W=gap(k)
        mu,b,tag=profile(k)
        directional_speed=abs(2*(low+W/2)*g/3)
        # Illustrative local crossing in a background with Omega/T
        # proportional to n^{1/3}, treating f and band shape as frozen.
        rows.append(dict(assumed_coherent_energy_fraction=f,background_mass_over_T=bg,Omega_over_T=O,k2=k,
                         relative_band_width=float(W/(low+W/2)),
                         local_peak_rate_in_Href=float(O*mu/H),
                         assumed_band_sweep_per_log_a=float(directional_speed),
                         estimated_amplitude_log_gain_in_one_crossing=float((O/H)*b/directional_speed),
                         method=tag))
    mu,b,_=profile(.002)
    W=gap(.002)[1]
    assert abs(mu/(W/8)-1)<.008
    assert abs(b/(np.pi*W**2/32)-1)<.008
    return dict(x=x,mass_ratio=s['mass_ratio'],source_dlogN_dloga=g,rows=rows,
                checks=dict(small_k2_numeric_to_analytic_peak_ratio=mu/(W/8),
                            small_k2_numeric_to_analytic_area_ratio=b/(np.pi*W**2/32)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    d={}
    for x in a.spectra:
        path=Path(x)
        raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
        d[path.name]=one(json.loads(raw))
    print(json.dumps(d,indent=2,allow_nan=False))
