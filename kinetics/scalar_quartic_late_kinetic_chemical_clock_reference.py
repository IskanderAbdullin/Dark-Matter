"""Separate vacuum quartic elastic encounters and upper chemical rates on held shapes."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre
from archived_collisionless_streaming_reference import Background
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import mp
from warm_prompt_reheating_expanding_slow_reference import visible
from scalar_quartic_vacuum_multiplication_bound_reference import global_upper


def rate(T,phase,radial,angular,mass):
    x,w=roots_legendre(radial)
    e=np.r_[0.,np.geomspace(2e-6,70.,17)]
    dx=(e[1:]-e[:-1])[:,None]/2
    p=(T*((e[1:]+e[:-1])[:,None]/2+dx*x)).ravel()
    E=np.hypot(p,mass)
    fp=phase*np.exp(-p*p/(E+mass)/T)
    q=(T*dx*w).ravel()*p*p*fp/(1.-fp)/(2*np.pi**2)
    uu,v=roots_legendre(angular)
    avg=np.zeros((len(p),len(p)))
    pp=2*p[:,None]*p[None,:]
    ep=2*E[:,None]*E[None,:]+2*mass*mass
    if mass:
        for u,a in zip(uu,v,strict=True):
            s=ep-u*pp
            avg+=a/2*np.sqrt(np.maximum(1.-4*mass*mass/s,0.))
    else:
        avg+=1.
    frac=(q/E)@avg@(q/E)/q.sum()**2
    return q.sum(),frac/(64*np.pi)


def run(eos,archive):
    bg=Background(eos)
    with gzip.open(archive,'rt') as fi:
        path=json.load(fi)['arms'][-1]['trajectory']
    couplings=(7.65175e-8,1.08515714e-6)
    rows=[]
    for st in (path[0],path[len(path)//2],path[-1]):
        T=st['T_visible']
        Td=st['light_scalar_T'][1];fug=st['phase_factors'][1]
        H=np.pi/np.sqrt(45)*T*T/MPL
        _,sv=visible(bg,T)
        n,low=rate(Td,fug,10,14,mp)
        n2,hi=rate(Td,fug,17,24,mp)
        assert abs(n/n2-1)<1e-7 and abs(hi/low-1)<2e-5
        assert abs(n2/(st['yields'][1]*sv)-1)<2e-7
        F=fug/(1-fug)
        rev=np.exp(2*(mp/Td+np.log(fug)))
        assert rev<1.
        # Exact massless MB is known analytically; use an extremely dilute BE.
        _,mb=rate(Td,1e-7,14,16,0.)
        assert abs(mb*(256*np.pi*Td*Td)-1)<1e-6
        ans=[]
        for L in couplings:
            # Photon-only H floor does NOT give a lower bound on the
            # rate/H. Store normalized diagnostic rates with clear labels.
            chem=n2*global_upper(L)*(1+F)**4/H
            kin=n2*hi*L*L/H
            ans.append(dict(lambda_quartic=L,elastic_unstimulated_rate_over_Hph=kin,
                            forward_net_count_gain_upper_over_Hph=chem,
                            reverse_net_count_loss_upper_over_Hph=chem*rev,
                            chemical_over_unstimulated_elastic_rate_upper=chem/kin))
        rows.append(dict(visible_T_GeV=T,scalar_T_GeV=Td,
                         reverse_to_forward_detailed_balance_for_held_BE=rev,
                         alternative_lambdas=ans))
    return dict(scope="Instantaneously held Bose shapes, tree-level real-scalar quartic only. "
                      "Elastic gross rate without final stimulation, not a relaxation eigenvalue; "
                      "Hph is a photon-only clock, not measured H. No change to the original evolution "
                      "or assertion that scalar scattering establishes equilibrium.",
                sampled=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','archive'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.archive),indent=2,allow_nan=False))
