"""Compare equal-count Bose spectra at different chemical-equilibrium mass ratios."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import zeta, roots_legendre
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import mp
from warm_prompt_reheating_expanding_slow_reference import visible


def I(power,r,upper=130.):
    return quad(lambda u:u**power/np.expm1(np.hypot(u,r)),0.,upper,
                epsabs=5e-15,epsrel=2e-11)[0]


def run(eos,history):
    bg=Background(eos)
    with gzip.open(history,'rt') as f:
        st=json.load(f)['arms'][-1]['trajectory'][-1]
    T=st['T_visible']
    a=bg.scale(T)
    _,s=visible(bg,T)
    conv=conversion(eos)['Omega_h2_per_GeV_times_yield']
    Y=.06/(conv*mp)
    q0=bg.T0
    nc=Y*s*a**3
    physical_n=nc/a**3
    theta=lambda r:(2*np.pi**2*nc/I(2,r))**(1/3)
    assert abs(I(2,0)/(2*zeta(3))-1)<2e-12
    rad=[]
    for r in (0.,.01,.1,1.,3.,6.,10.):
        th=theta(r)
        series=dict(m_over_Tfreeze=r,comoving_T_over_Tgamma0=th/q0,
                    mean_q_over_Tgamma0=th/q0*I(3,r)/I(2,r),
                    rms_q_over_Tgamma0=th/q0*np.sqrt(I(4,r)/I(2,r)),
                    fraction_q_over_Tgamma0_below_point1=I(2,r,.1*q0/th)/I(2,r),
                    f_at_rest=None if r==0 else 1/np.expm1(r))
        rad.append(series)
    for x,y in zip(rad[:-1],rad[1:],strict=True):
        r1,r2=x['m_over_Tfreeze'],y['m_over_Tfreeze']
        t1,t2=theta(r1),theta(r2)
        qc=np.sqrt((r2*r2-r1*r1)/(t1**-2-t2**-2))
        # The difference in squared exponential arguments changes sign once.
        assert np.hypot(.51*qc/t1,r1)<np.hypot(.51*qc/t2,r2)
        assert np.hypot(1.7*qc/t1,r1)>np.hypot(1.7*qc/t2,r2)
        assert y['fraction_q_over_Tgamma0_below_point1']<x['fraction_q_over_Tgamma0_below_point1']
        y['one_crossing_q_over_Tgamma0_relative_to_previous_row']=qc/q0
    r=brentq(lambda z:theta(z)*z-a*mp,1.e-9,1.,xtol=2e-15)
    Td=mp/r
    uu,ww=roots_legendre(24)
    edge=np.r_[0.,np.geomspace(1e-6,84.,31)]
    lo,hi=edge[:-1,None],edge[1:,None]
    x=((lo+hi)/2+(hi-lo)/2*uu).ravel()
    weights=((hi-lo)/2*ww).ravel()*x*x/np.expm1(np.hypot(x,r))
    rr=weights*np.hypot(x,r)
    e_required=Td**4*rr.sum()/(2*np.pi**2)
    assert abs(Td**3*weights.sum()/(2*np.pi**2)/physical_n-1)<5e-9
    stream=(weights@bg.streams(x*Td*a,mp,a,bg.aeq,8))/weights.sum()
    sd=sum(st['energy_per_visible_entropy'])*s
    def energy(temp):
        r=mp/temp
        return temp**4/(2*np.pi**2)*quad(
            lambda x:x*x*np.hypot(x,r)/np.expm1(np.hypot(x,r)),
            0.,130.,epsabs=1e-13)[0]
    T_all=brentq(lambda q:energy(q)/sd-1,Td/5,Td*6)
    n_all=T_all**3*I(2,mp/T_all)/(2*np.pi**2)
    return dict(scope="One normal-phase real scalar, ideal Bose chemical equilibrium with mu=0 "
                      "at a hypothetical sudden freeze of the comoving count, no condensate, subsequent "
                      "momentum redshift, no changed visible entropy, comparison fixed mass and count. "
                      "No rate calculation or guarantee of early equilibration or subsequent decoupling.",
                target_Y=Y,dimensionless_comoving_shape_ledger=rad,
                modeled_anchor_counterfactual=dict(visible_T_GeV=T,
                 required_dark_temperature_GeV=Td,
                 required_scalar_energy_share_of_recorded_dark_total=e_required/sd,
                 required_scalar_energy_over_recorded_scalar_energy=e_required/(s*st['energy_per_visible_entropy'][1]),
                 average_path_from_anchor_in_fixed_background_Mpc=float(stream),
                 all_recorded_dark_energy_to_scalar_equilibrium_T_GeV=T_all,
                 all_energy_thermal_Y=n_all/s,
                 all_energy_thermal_scalar_omega_if_counts_freeze=n_all/s*conv*mp))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('history',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.history),indent=2,allow_nan=False))
