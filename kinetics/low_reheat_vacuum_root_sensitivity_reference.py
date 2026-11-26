"""Conditional perturbations of the low-reheat vacuum number normalization."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from empty_lepton_pair_spectrum_reference import Source
from late_radiation_background_reference import TabulatedRadiation
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M,m,mp


class NumberFold:
    def __init__(self,table,masses,order,time):
        self.bg=TabulatedRadiation(table,.009)
        self.mass=np.asarray(masses)
        self.src=[Source(order,dark_mass=x) for x in (*masses,M)]
        self.x,self.w=roots_legendre(time)
        lo,hi=np.log([.001,.025])
        self.edges=np.unique(np.clip(np.r_[lo,hi,self.bg.sg.x],lo,hi))
        panel=[self.panel(a,b) for a,b in zip(self.edges[:-1],self.edges[1:])]
        self.cum=np.cumsum(np.vstack([np.zeros(len(self.src)),panel]),axis=0)

    def rate_vec(self,u):
        T=np.exp(u)
        s=2*np.pi**2/45*np.exp(self.bg.entropy_log(u))
        H=np.sqrt(np.pi**2/90*np.exp(4*u+self.bg.rg(u)))/MPL
        measure=(3.+self.bg.sg.derivative()(u))/(3*H*s)
        return 2*measure*np.array([x.event_moments(T)[3].sum() for x in self.src])

    def panel(self,a,b):
        if b==a:
            return np.zeros(len(self.src))
        u=(a+b)/2+(b-a)/2*self.x
        return (self.w@np.asarray([self.rate_vec(x) for x in u]))*(b-a)/2

    def at(self,T):
        u=np.log(T)
        j=np.searchsorted(self.edges,u,side="right")-1
        assert 0<=j<len(self.edges)-1
        return self.cum[j]+self.panel(self.edges[j],u)


def run(eos):
    masses=m*np.array([.9,1.0,1.1])
    nfactor=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    target=.120
    controls=[]
    for order,t in (((6,15,12),4),((10,25,18),6)):
        K=NumberFold(eos,masses,order,t)
        rows=[]
        for j,mj in enumerate(masses):
            def f(T,fac=1.):
                R=K.at(T)
                return fac*fac*nfactor*np.array([mp*R[-1],mj*(R[j]+R[-1])])
            def root(fac):
                return brentq(lambda T:f(T,fac)[1]-target,.004,.016,xtol=1e-12)
            nom=root(1.)
            g=K.at(nom)
            gp=K.rate_vec(np.log(nom))
            log_slope=(gp[j]+gp[-1])/(g[j]+g[-1])
            h=1e-4
            dif=(np.log(f(nom*np.exp(h))[1])-np.log(f(nom*np.exp(-h))[1]))/(2*h)
            assert abs(log_slope/dif-1)<2e-7
            rminus,rplus=root(1.-h),root(1.+h)
            epsd=(np.log(rplus)-np.log(rminus))/(np.log(1.+h)-np.log(1.-h))
            assert abs(epsd/(-2/log_slope)-1)<2e-7
            rows.append(dict(stable_mass_GeV=mj,root_visible_reheat_GeV=nom,
                scalar_direct_prompt_omega_h2_at_root=f(nom)[0],
                local_log_count_slope_to_log_TR=log_slope,
                local_log_root_slope_to_log_portal_mixing=-2/log_slope,
                finite_difference_TR_slope=dif,finite_difference_mixing_slope=epsd,
                stable_omega_at_minus_plus1pct_TR=[f(nom*x)[1] for x in (.99,1.01)],
                roots_for_mixing_multipliers_0p9_1p1_GeV=[root(.9),root(1.1)]))
        controls.append(dict(order=order,time_order=t,rows=rows))
    num=max(abs(a["root_visible_reheat_GeV"]/b["root_visible_reheat_GeV"]-1.)
       for a,b in zip(controls[0]["rows"],controls[1]["rows"]))
    assert num<5e-5
    return dict(reference_stable_Omega_h2=target,eos_sha256=hashlib.sha256(Path(eos).read_bytes()).hexdigest(),
                max_root_refinement_difference=num,controls=controls,
                scope="Vacuum final states, exact FD incoming leptons, prompt heavy decay counts, prescribed thermal radiation from instant reheating to 1 MeV. Mixing variation multiplies the fixed spectral normalization by epsilon squared at fixed widths and masses. No feedback, inverse depletion, dark reactions, dominant scalar or fit.")


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("eos"); a=p.parse_args()
    print(json.dumps(run(a.eos),indent=2,allow_nan=False))
