"""Future-timelike denominator floors for on-shell charged external lines."""
import gzip
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_frozen_longitudinal_pole_reference import L
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF, Species
from finite_mass_qed_collective_pole_pair_threshold_reference import run as bounds


def denominators(m,p,mu,k,c):
    E=np.hypot(m,p)
    W=np.hypot(mu,k)
    lower=(m*m*W*W+mu*mu*p*p)/(E*W+p*k)
    stable=lower+p*k*(1-c)
    return mu*mu+2*stable,mu*mu-2*stable,stable


def one(g,d):
    frozen=bounds(g,d)
    low=frozen['plasma_gap']
    high=frozen['transverse_asymptotic_invariant_mass']
    p=d['x_end']
    species=[Species(m,d[f'endpoint_{tag}_momentum_edges'],d[f'endpoint_{tag}_occupations'])
             for m,tag in ((p,'parent'),(p*d['mass_ratio'],'fermion'))]
    floors=[]; equality=[]; direct=[]; Lnear=[]
    if not high<min(c.m for c in species):
        raise ArithmeticError('This archiving bound assumes monotone magnitude polynomials.')
    for S in species:
        m=S.m
        floors.append(dict(mass=m,plus_uniform_floor=low*(2*m+low),
                           minus_magnitude_uniform_floor=low*(2*m-low)))
        for beta in (0.,.01,.27,.7,.92,.98,.999,.99999):
            T=np.sqrt(sum(c.field(beta)[0] for c in species))
            k=T*beta/np.sqrt((1-beta)*(1+beta))
            if beta==0.:
                C=low
            else:
                C=np.sqrt(PREF*sum(c.W@L(c.v*beta)[0] for c in species)
                          *((1-beta)*(1+beta))/(beta*beta))
            kL=C*beta/np.sqrt((1-beta)*(1+beta))
            for mu,q,name in ((T,k,'T'),(C,kL,'L')):
                p0=m*q/mu
                a,b,x=denominators(m,p0,mu,q,1.)
                equality.append(abs(x/(m*mu)-1.))
                if name=='L' and beta>=.98:
                    Lnear.append(dict(m=m,beta=beta,mu=mu,comoving_minus_magnitude=-b))
                for mom in (0.,.5*p0,p0,2*p0,.01,.25,1.,10.,50.):
                    for c in (-1.,-.3,0.,.4,.99,1.):
                        plus,minus,X=denominators(m,mom,mu,q,c)
                        if plus<mu*(2*m+mu)*(1-4e-13) or -minus<mu*(2*m-mu)*(1-4e-13):
                            raise ArithmeticError('Analytic future-timelike bound failed')
                        if name=='T' and (plus<floors[-1]['plus_uniform_floor']*(1-4e-13)
                                          or -minus<floors[-1]['minus_magnitude_uniform_floor']*(1-4e-13)):
                            raise ArithmeticError('An archive transverse denominator lost its uniform margin')
                        via=np.hypot(m,mom)*np.hypot(mu,q)-mom*q*c
                        direct.append(abs(X-via)/(m*mu))
    if max(equality)>2e-13 or max(direct)>2e-8:
        raise ArithmeticError((max(equality),max(direct)))
    return dict(grid=g,transverse_floors=floors,comoving_sharpness_max_rel=max(equality),
                direct_dot_absolute_defect_in_m_mu_units=max(direct),direct_checks=len(direct),
                longitudinal_light_cone_examples=Lnear)


if __name__=='__main__':
    p=Path(__file__).parent/'results'
    out=[]
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:
            out.append(one(g,json.load(f)))
    print(json.dumps(dict(scope='Timelike pole kinematics and bare denominators only; no cut matching or rate.',data=out),
                     indent=2,allow_nan=False))
