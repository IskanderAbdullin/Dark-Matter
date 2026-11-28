"""Coldest ideal Bose normal phase at fixed count, including a NR limiting control."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gamma, zeta
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import mp
from warm_prompt_reheating_expanding_slow_reference import visible


def fweight(u,r,power):
    arg=u*u/(np.hypot(u,r)+r) if r else u
    return u**power*(np.exp(-arg) if arg>80 else 1/np.expm1(arg))


def I(r,power,cap=None):
    scale=np.sqrt(max(1.,2*r))
    if r>0 and cap is None:
        # In the nonrelativistic regime use u~sqrt(2r), including the
        # crossover and tails rather than relying on a fixed cutoff.
        lim=200+40*scale
    else:
        lim=180.
    if cap is not None:
        lim=cap
    return quad(lambda u:fweight(u,r,power),0.,lim,epsabs=1e-12,epsrel=5e-11)[0]


def J(power,cap=None):
    if cap is None:
        return gamma((power+1)/2)*zeta((power+1)/2)/2
    return quad(lambda v:v**power/np.expm1(v*v),0.,min(cap,24.),epsabs=2e-13)[0]


def run(eos,history):
    bg=Background(eos)
    with gzip.open(history,'rt') as f:
        st=json.load(f)['arms'][-1]['trajectory'][-1]
    T=st['T_visible']
    a=bg.scale(T)
    _,s=visible(bg,T)
    Y=.06/(mp*conversion(eos)['Omega_h2_per_GeV_times_yield'])
    qg=bg.T0
    nc=Y*s*a**3
    theta=lambda r:(2*np.pi**2*nc/I(r,2))**(1/3)
    theta_nr=(2*np.pi**2*nc/J(2))**(1/3)
    nr=dict(mean_q_over_Tgamma0=theta_nr/qg*J(3)/J(2),
            rms_q_over_Tgamma0=theta_nr/qg*np.sqrt(J(4)/J(2)),
            cold_frac_point1=J(2,.1*qg/theta_nr)/J(2),
            cold_frac_point01=J(2,.01*qg/theta_nr)/J(2))
    R=brentq(lambda r:r*theta(r)-a*mp,1.e-9,1.,xtol=2e-15)
    data=[]
    for r in (0.,R,.1,1.,3.,10.,100.,10000.):
        t=theta(r)
        it=I(r,2)
        row=dict(m_over_T=r,scale_q_over_Tgamma0=t/qg,
                 rms_q_over_Tgamma0=t/qg*np.sqrt(I(r,4)/it),
                 mean_q_over_Tgamma0=t/qg*I(r,3)/it,
                 cold_frac_point1=I(r,2,.1*qg/t)/it,
                 cold_frac_point01=I(r,2,.01*qg/t)/it)
        if data:
            v=data[-1]
            assert row['mean_q_over_Tgamma0']<v['mean_q_over_Tgamma0']
            assert row['cold_frac_point1']>v['cold_frac_point1']
        data.append(row)
    assert abs(data[-1]['cold_frac_point1']/nr['cold_frac_point1']-1)<2e-4
    assert abs(data[-1]['rms_q_over_Tgamma0']/nr['rms_q_over_Tgamma0']-1)<2e-4
    flags=[]
    for k in (.01,.1):
        for f in (.5,.9):
            def all_cold(c):
                n=1-c
                return c+n*J(2,k*qg/(theta_nr*n**(1/3)))/J(2)
            c=brentq(lambda c:all_cold(c)-f,0.,1.-1e-9)
            flags.append(dict(q_cut_over_Tgamma0=k,target_cold_fraction=f,
                              minimum_condensate_fraction_if_excited_component_is_normal_thermal=c))
    return dict(scope="Ideal one-species Bose equilibrium, including limiting excited state at mu=m, "
                      "fixed comoving total scalar count at requested final rest abundance, no "
                      "effective mass shift. A hypothetical condensate treated as zero momentum only "
                      "for geometric counts; neither dynamics nor free-streaming from a physical freeze.",
                target_Y=Y,anchors=data,ideal_NR_saturated_lower_warmth=nr,
                cold_population_condensate_examples=flags)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('history',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.history),indent=2,allow_nan=False))
