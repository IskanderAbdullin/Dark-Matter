"""Hard central-angle Abelian photon pair-production phase-space reference."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit, roots_legendre


def central(s,m,g,amax):
    r=np.minimum(1.,4*m*m/s)
    t=amax*np.sqrt(1-r)
    F=(4+4*r-2*r*r)*np.arctanh(t)-2*t-2*r*r*t/(1-t*t)
    return g**4*F/(8*np.pi*s)


def control():
    m=.34; s=2.
    r=4*m*m/s
    beta=np.sqrt(1-r)
    for cut in (.6,.8,.95,1.):
        v=quad(lambda t:beta*(4*(1+r)/(1-(beta*t)**2)-2-
                              4*r*r/(1-(beta*t)**2)**2)/(16*np.pi*s),
               -cut,cut,epsabs=1e-13)[0]
        assert abs(v/central(s,m,1.,cut)-1)<2e-12
    bw=((3-beta**4)*np.log((1+beta)/(1-beta))-2*beta*(2-beta*beta))/(8*np.pi*s)
    assert abs(central(s,m,1.,1.)/bw-1)<2e-12


def calc(m,g,cut,order):
    z,w=roots_legendre(order)
    p=.5+(z+1)*19.75
    pop=2*p*p/np.expm1(p)*w*19.75/(2*np.pi**2)
    s0=4*p[:,None]*p[None,:]
    a=np.minimum(1.,np.maximum(1.,4*m*m)/s0)
    out=np.zeros_like(s0)
    for x,y in zip(z,w,strict=True):
        t=a+.5*(1-a)*(x+1)
        out+=y*(1-a)*t*central(s0*t,m,g,cut)
    # 1/2 for identical incident photons, multiplied by 2 particles.
    return pop@out@pop


def one(r):
    T=r['joint_temperature_over_T']
    M=r['mass_ratio']*r['x']/T
    mu=r['daughter_mu_over_T']/T
    N=r['initial_daughter_number_over_T3']
    gate=r['thermal_daughter_number_at_zero_stimulation_over_T3']
    H=np.sqrt(np.pi**2*106.75/90)*(86./r['x'])/2.435e18
    data=[]
    for g in (.0003,.001,.003,.01):
        mass=np.hypot(M,g/2)
        rates=np.array([calc(mass,g,.8,i) for i in (40,80)])
        assert abs(rates[0]/rates[-1]-1)<2e-4
        gain=rates[-1]*T**4/H
        net=(1-np.exp(2*mu))*(1-expit(mu-mass))**2*gain
        data.append(dict(g=g,unblocked_hard_pair_particles_per_log_a_over_T3=gain,
                         net_lower_limited_vacuum_initial_slice=net,
                         initial_target_daughter_number=N,number_gap_to_fixed_population_gate=gate-N))
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('gates')
    a=p.parse_args()
    d=json.loads(Path(a.gates).read_text())
    control()
    print(json.dumps({k:one(v) for k,v in d.items()},indent=2,allow_nan=False))
