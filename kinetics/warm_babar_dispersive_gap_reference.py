"""Bound the 0.30–sqrt(0.10) GeV pion gap from a positive BABAR moment."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre, kve

from late_radiation_background_reference import MPL, TabulatedRadiation
from measured_partial_hadron_reheat_fold_reference import integrands
from three_species_decay_onsager_reference import m, M
from timelike_portal_onsager_coverage_reference import V, alpha, factor, mpi, width


POINT=4*np.pi*alpha**2/3  # sigma_point = POINT/s, sigma in GeV^-2.
MU=.1056583755


def mu_kernel_lower(s, n=1000):
    e=np.linspace(0.,1.,n+1)
    lo,hi=e[:-1],e[1:]
    # Lower bound numerator and upper bound denominator by endpoints;
    # safe also on either end of [0,1], when the lower numerator is 0.
    return (np.diff(e)*(lo**2)*(1-hi)) @ (
        1/(hi[None,:]**2+(1-lo[None,:])*np.asarray(s)[:,None]/MU**2)).T


def thermal_window_upper(eos,n_s,n_t):
    l,h=np.linspace(.30**2,.10,n_s+1)[:-1],np.linspace(.30**2,.10,n_s+1)[1:]
    El,Eh=np.sqrt(l),np.sqrt(h)
    H_kernel=mu_kernel_lower(h)/(4*np.pi**3)
    assert min(H_kernel)>0
    x,w=roots_legendre(120)
    x=(x+1)/2
    exact=((x*x*(1-x))*w/2) @ (1/(x*x+(1-x)*h[:,None]/MU**2)).T/(4*np.pi**3)
    assert min(exact/H_kernel)>1. and max(exact/H_kernel)<1.02
    A=factor*h**2.5/(POINT*((h-V**2)**2+(V*width)**2))
    A*=sum((h+2*d*d)*np.sqrt(1-4*d*d/h) for d in (m,M))
    T_check=.05
    zz=Eh/T_check
    exact_local=A*np.exp(-zz)*(2*T_check**2*kve(1,zz)+T_check*Eh*kve(0,zz))
    norm=np.asarray([integrands(T_check,[(np.array([s]),np.array([s/POINT]),np.zeros(1))])[0][:,1].sum()
                     for s in h])
    check=max(abs(exact_local/norm-1.))
    assert check<5e-14

    T0,b= .05,500.
    R=TabulatedRadiation(eos,T0)
    Hph=np.pi*T0**2/(np.sqrt(45)*MPL)
    ll,hh=np.log(R.state(b)[0]),R.initial
    p=R.sg.x[(R.sg.x>ll)&(R.sg.x<hh)]
    u=np.unique(np.r_[np.linspace(ll,hh,n_t+1),p])[::-1]
    gs,gr=R.sg(u),R.rg(u)
    t=(R.s0-3*u-gs)/3
    hr=np.exp(.5*(np.minimum(gr[:-1],gr[1:])-np.log(2))+
              (2/3)*(gs[0]-np.maximum(gs[:-1],gs[1:]))-1e-12)
    assert np.all(np.diff(t)>0.) and abs(t[-1]-np.log(b))<1e-12
    T=np.exp(u[:-1,None])
    BE=(-1/np.expm1(-mpi/T))**2
    z=El[None,:]/T
    K0,K1=[kve(n,z)*np.exp(-z) for n in (0,1)]
    Q=BE*(2*T*T*K1+T*Eh[None,:]*K0)
    J=np.diff(np.exp(5*t))/(5*Hph*hr) @ Q
    # For z>=100, K_nu(z) <= exp(-100) sqrt(pi/198)*(100/z)^8.
    smax=np.exp((gs[0]-np.log(2))/3+1e-12)*T0
    assert min(El)/(smax/b)>100.
    Tcold=smax/b
    f=(-1/np.expm1(-mpi/Tcold))**2*np.exp(-100)*np.sqrt(np.pi/198)/Hph
    J_cold=f*(100/El)**8*(2*smax**10/(5*b**5)+Eh*smax**9/(4*b**4))
    gain=A*(J+J_cold)/H_kernel
    j=np.argmax(gain)
    return dict(s_cells=n_s,temperature_cells=len(t)-1,
                Q_over_positive_moment_upper_GeV4=float(gain[j]*(1+1e-11)),
                controlling_s_interval_GeV2=[float(l[j]),float(h[j])],
                max_gauss_over_kernel_lower=float(max(exact/H_kernel)),
                normalization_vs_previous_hadronic_kernel_max_defect=float(check),
                max_cold_share=float(max(J_cold/(J+J_cold))))


def run(eos,previous):
    old=json.loads(Path(previous).read_bytes())
    pion=next(z["comparison"] for z in old["cases"]
              if z["policy"]=="same_band_with_massive_initial_Bose_upper")
    budget=pion["extra_energy_after_including_pions_to_touch_half_GeV4"]
    base=pion["phased_measured_pion_and_lepton_scalar_Omega_ceiling"]
    folds=[thermal_window_upper(eos,n,t) for n,t in ((40,550),(90,1200))]
    # BABAR Table I. These are conditional input choices, not frequentist
    # confidence limits, and the two measurements are not combined here.
    variants=[]
    for label,value,stat,syst in (("BABAR 2026",60.7,1.7,1.4),("BABAR combination",58.2,.86,.78)):
        for multiples in (1,2):
            proposed=(value+multiples*(stat+syst))*1e-10
            gain=folds[-1]["Q_over_positive_moment_upper_GeV4"]*proposed
            share=gain/budget
            assert 0.<share<1.
            # The positive quadratic in earlier deposited energy is convex.
            # Its secant is an independent conservative upper at the gap budget.
            secant=base+share*(.060-base)
            variants.append(dict(data=label,linear_error_multiple=multiples,
                                 input_whole_dispersion_moment=proposed,
                                 gap_energy_upper_GeV4=gain,share_of_available_extra_energy=share,
                                 remaining_budget_GeV4=budget-gain,
                                 scalar_Omega_upper_with_this_gap_if_no_else=secant))
    return dict(scope="Positive moment bound for nonoverlap 0.30–sqrt(.10)GeV, effective "
                      "massive two-charged-pion vacuum kernel with bounded Bose weights. "
                      "The experimental moment cap is a conditional input. No data assumed "
                      "under 0.30 GeV or other hadronic sources resolved.",
                source="https://arxiv.org/html/2611.04596v1#S0.T1",gap_s_GeV2=[.09,.10],
                prior_remaining_after_measured_bosonic_band_GeV4=budget,
                prior_scalar_Omega_upper=base,folds=folds,
                moment_variants=variants)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("prior")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior),indent=2,allow_nan=False))
