"""All-branch hard-loop invariant bounds versus the exact free pair threshold."""
import gzip
import json
from pathlib import Path

import numpy as np

from finite_mass_qed_frozen_longitudinal_pole_reference import L
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import PREF, Species, gl


def run(grid,d):
    m=d['x_end']
    args=[(M,d[f'endpoint_{tag}_momentum_edges'],d[f'endpoint_{tag}_occupations'])
          for M,tag in ((m,'parent'),(m*d['mass_ratio'],'fermion'))]
    species=[Species(*x) for x in args]
    mass_threshold=2*min(c.m for c in species)
    x,w=gl(56)
    transverse_ibp=0.
    longitudinal_ibp=np.zeros(26)
    for c in species:
        a,b=c.e[:-1,None],c.e[1:,None]
        p=(a+b)/2+(b-a)*x/2
        E=np.hypot(c.m,p)
        v2=(p/E)**2
        dw=PREF*c.f*(b-a)[:,0]/2
        transverse_ibp+=dw@((p*p/E)@w)
        for j,n in enumerate(range(1,27)):
            dA=E*v2**n*(1-((2*n-1)/(2*n+1))*v2)
            longitudinal_ibp[j]+=dw@(dA@w)
    pi_inf=sum(c.field(1.)[0] for c in species)
    wp2=sum(c.G for c in species)
    defect_T=abs(pi_inf/transverse_ibp-1.)
    defect_L=abs(wp2/longitudinal_ibp[0]-1.)
    if max(defect_T,defect_L)>2e-12 or np.max(np.diff(longitudinal_ibp))>=0.:
        raise ArithmeticError((grid,defect_T,defect_L,longitudinal_ibp.tolist()))
    betas=np.r_[0.,np.geomspace(1e-5,.8,45),np.linspace(.81,1.,50)]
    T,Lv=[],[]
    for beta in betas:
        T.append(sum(c.field(beta)[0] for c in species))
        Lv.append(wp2 if beta==0. else PREF*sum(c.W@L(c.v*beta)[0] for c in species)
                  *((1-beta)*(1+beta))/(beta*beta))
    T,Lv=np.array(T),np.array(Lv)
    if np.any(np.diff(T)<-1e-17) or np.any(np.diff(Lv)>1e-17) or abs(Lv[-1])>1e-17:
        raise ArithmeticError('Invariant-mass monotonicity sweep failed')
    if max(pi_inf,wp2)>=mass_threshold**2:
        raise ArithmeticError('A threshold can open in the held hard loop')
    return dict(grid=grid,free_pair_minimum_mass=mass_threshold,plasma_gap=np.sqrt(wp2),
                transverse_asymptotic_invariant_mass=np.sqrt(pi_inf),
                longitudinal_max_pair_threshold_ratio=np.sqrt(wp2)/mass_threshold,
                transverse_max_pair_threshold_ratio=np.sqrt(pi_inf)/mass_threshold,
                transverse_asymptotic_invariant_ibp_relative_defect=defect_T,
                longitudinal_plasma_gap_ibp_relative_defect=defect_L,
                longitudinal_ibp_coefficients=longitudinal_ibp.tolist(),
                beta_scan_size=len(betas),
                smallest_transverse_invariant_mass=np.sqrt(T.min()))


if __name__=='__main__':
    r=Path(__file__).parent/'results'
    data=[]
    for g in ('fine','third'):
        with gzip.open(r/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as o:
            data.append(run(g,json.load(o)))
    print(json.dumps(dict(scope='Frozen timelike collisionless hard-loop pole bound vs exact free charged-pair threshold.',
                          data=data),indent=2,allow_nan=False))
