"""Riemann upper source and count bound for the fixed-g instantaneous history."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import gamma, gammaincc, kve

from late_radiation_background_reference import TabulatedRadiation, MPL
from three_species_decay_onsager_reference import M, m, mp
from warm_late_visible_recycling_envelope_reference import (
    COUP2, LEPTONS, VECTOR_GAMMA, VECTOR_MASS, energy_rate_coefficient,
    temperature_time_upper,
)


GSTAR=10.75
MPLANCK=1.2209e19
S0=2891.2
RHOC_H2=1.05375e-5


def invariant_tail_integral(u):
    assert u>1.
    # Bound the Bessel functions and then bound each term by a gamma tail.
    return (np.sqrt(np.pi/(2*(u-1))) *
            (2*gamma(7)*gammaincc(7,u)+gamma(8)*gammaincc(8,u)) /
            (VECTOR_MASS*VECTOR_GAMMA)**2)


def upper_pair_Q(T,du,split=2.):
    terms=0.
    direct=0.
    for d in (m,M):
        for l in LEPTONS:
            lower = 2*max(d,l)/T
            upper = min(split/T, lower+55.)
            u0 = max(lower, upper)
            tail = invariant_tail_integral(u0)
            direct_tail = (np.sqrt(np.pi/(2*(u0-1))) * gamma(7) *
                           gammaincc(7,u0)/(VECTOR_MASS*VECTOR_GAMMA)**2)
            final_particles = 2 if d == m else 4
            direct += final_particles*direct_tail
            if upper <= lower:
                terms += tail
                continue
            z = np.linspace(lower,upper,int(np.ceil((upper-lower)/du))+1)
            sl=(T*z[1:])**2
            def phase(r):
                return (1+2*r)*np.sqrt(np.maximum(0.,1-4*r))
            # The thermal K functions are decreasing. All u powers,
            # phase factors and the sub-pole propagator increase.
            K0,K1 = [kve(n,z[:-1])*np.exp(-z[:-1]) for n in (0,1)]
            norm=1/((sl-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2)
            positive = z[1:]**6*phase(d*d/sl)*phase(l*l/sl)*norm*np.diff(z)
            terms += positive @ (2*K1+z[1:]*K0) + tail
            direct += final_particles*(positive @ K1)
    C=1.00000001*COUP2*T**8/(48*np.pi**5)
    return np.array([C*T*terms,C*direct])


def run_case(TR,K,ntime,du,stop=5.):
    H=1.66*np.sqrt(GSTAR)*TR*TR/MPLANCK
    NORM=mp*S0/(RHOC_H2*(2*np.pi**2/45*GSTAR*TR**3))
    t=np.linspace(0,stop,ntime+1)
    rates=np.array([upper_pair_Q(T,du) for T in TR*np.exp(-t[:-1])])
    assert np.all(np.diff(rates,axis=0)<=0.)
    factor=np.diff(np.exp(5*t))/(5*H)
    Qs=rates[:,0]*factor
    cumulative=np.cumsum(Qs)
    # Use a uniform massless tail after stop, letting the finite source
    # upper envelope remain available at arbitrarily large expansion.
    tail_rate=energy_rate_coefficient(.01045)[0]*TR**9
    q_tail=tail_rate*np.exp(-4*stop)/(4*H)
    total_energy=cumulative[-1]+q_tail
    heat_exp=288*COUP2*TR**5/(3*np.pi**3*VECTOR_MASS**4*H)
    recycle=(K/(M*M*H)*np.exp(2*heat_exp) *
             (np.diff(-np.exp(-t)) @ cumulative**2 + np.exp(-stop)*total_energy**2))
    direct=rates[:,1]@factor + q_tail/m
    assert M>=2*m and direct<total_energy/m
    return dict(TR_GeV=TR,time_panels=ntime,invariant_panel_width=du,
                conservative_energy_input_GeV4=float(total_energy),
                terminal_source_energy_tail_relative=float(q_tail/total_energy),
                heating_exponent=float(heat_exp),
                Omega_phi_direct_count_ceiling=float(direct*NORM),
                Omega_phi_recycle_ceiling=float(recycle*NORM),
                Omega_phi_total_ceiling=float((direct+recycle)*NORM))


def with_eos(TR,K,ntime,du,eos,stop=500.):
    rad=TabulatedRadiation(eos,TR)
    Tlo=rad.state(stop)[0]
    lr=rad.initial
    l0=np.log(Tlo)
    extra=rad.sg.x[(rad.sg.x>l0)&(rad.sg.x<lr)]
    u=np.unique(np.r_[np.linspace(l0,lr,ntime+1),extra])[::-1]
    lgs,lgr=rad.sg(u),rad.rg(u)
    t=(rad.s0-3*u-lgs)/3
    assert np.min(np.diff(t))>0 and abs(t[-1]-np.log(stop))<1e-12
    # Intervals never cross a PCHIP knot. Bound both EOS factors at
    # their independently worst endpoint instead of sampling H.
    Hph=np.pi/np.sqrt(45)*TR*TR/MPL
    floor=np.exp(.5*(np.minimum(lgr[:-1],lgr[1:])-np.log(2))
                 +(2/3)*(lgs[0]-np.maximum(lgs[:-1],lgs[1:]))-1e-12)
    rates=np.asarray([upper_pair_Q(T,du) for T in np.exp(u[:-1])])
    assert np.all(np.diff(rates,axis=0)<=0.)
    factor=np.diff(np.exp(5*t))/(5*Hph*floor)
    cumulative=np.cumsum(rates[:,0]*factor)
    over=np.exp((lgs[0]-np.log(2))/3+1e-12)
    late_pair=energy_rate_coefficient(.01045)[0]*(TR*over)**9/(4*Hph*stop**4)
    assert TR*over/stop < .01045
    etot=cumulative[-1]+late_pair
    heat=(288*COUP2/(np.pi**3*VECTOR_MASS**4))*temperature_time_upper(eos,TR,stop)[1]
    recyclable=K*np.exp(2*heat)/(M*M*Hph) * (
        (np.diff(-np.exp(-t))/floor) @ cumulative**2 + etot**2/stop)
    direct=rates[:,1]@factor+late_pair/m
    assert direct<etot/m
    entropy=2*np.pi**2/45 * np.exp(lgs[0])*TR**3
    norm=mp*S0/(RHOC_H2*entropy)
    return dict(TR_GeV=TR,time_panels=len(t)-1,invariant_panel_width=du,
                conservative_energy_input_GeV4=float(etot),
                terminal_source_energy_tail_relative=float(late_pair/etot),
                heating_exponent=float(heat),
                Omega_phi_direct_count_ceiling=float(direct*norm),
                Omega_phi_recycle_ceiling=float(recyclable*norm),
                Omega_phi_total_ceiling=float((direct+recyclable)*norm))


def run(kernel,eos=None):
    data=json.loads(Path(kernel).read_bytes())
    k=data["simultaneous_energy_density_coefficient_GeVminus2"]
    controls=[run_case(.05,k,n,d) for n,d in ((250,.03),(500,.015),(1000,.0075),(1500,.005))]
    assert controls[0]["Omega_phi_total_ceiling"]>controls[-1]["Omega_phi_total_ceiling"]
    eos_cases=[with_eos(.05,k,n,d,eos) for n,d in ((500,.015),(1200,.006))] if eos else []
    if eos_cases:
        assert eos_cases[-1]["Omega_phi_total_ceiling"]<eos_cases[0]["Omega_phi_total_ceiling"]
    return dict(scope="Initially empty fixed-RD g*=10.75 instantaneous vacuum "
                      "charged-lepton benchmark, plus an optional tabulated visible EOS clock. "
                      "No other channels or alternative reheating histories tested.",
                K_GeVminus2=k,cosmology=dict(gstar=GSTAR,Planck= MPLANCK,S0=S0,rhoc_over_h2=RHOC_H2),
                fixed_g_controls=controls,tabulated_EOS_controls=eos_cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("kernel")
    p.add_argument("--eos")
    a=p.parse_args()
    print(json.dumps(run(a.kernel,a.eos),indent=2,allow_nan=False))
