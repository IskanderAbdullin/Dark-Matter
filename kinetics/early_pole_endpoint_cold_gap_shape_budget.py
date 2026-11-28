"""Separate the cold scalar count and cold shape effects in the nodal Hartree gap."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def card(c,n):
    g=c['grid']
    cut=c['cut']
    if (g,c['scalar_div'],c['charge_div'],c['x_end'])!=('fine',n,2,.05):
        raise ArithmeticError('Use the matched fine S2 and S4 states')
    d=ScalarSubshellChart(load(g)[0],SPECS[g],2,n,cut)
    f=np.asarray(c['endpoint_occupation'])
    if c['counts']!=[d.ns,d.nf,d.np] or f.shape!=(d.n,) or not np.all(np.isfinite(f)) or min(f)<0:
        raise ArithmeticError('Invalid scalar shape')
    return d,f,d.state(c['x_end'],f)[0]


def evaluate(d,f,m,qcut):
    q=d.ds.p
    N=d.ds.w*f[:d.ns]
    cold=d.ds.ed[1:]<=qcut
    E=np.hypot(q,m)
    gap_term=N/E
    count=fsum(N[cold])
    remainder=fsum((N*q*q/(m*E*(E+m)))[cold])
    warm=fsum(gap_term[~cold])
    total=fsum(gap_term)
    if abs((count/m-remainder+warm)/total-1)>3e-14:
        raise ArithmeticError('Nodal gap ledger does not close')
    return count,remainder,warm,total


def run(states):
    a,b=[card(x,n) for x,n in zip(states,(2,4),strict=True)]
    qcut=states[0]['cut']
    if states[1]['cut']!=qcut or states[0]['steps']!=states[1]['steps']:
        raise ArithmeticError('Incompatible scalar geometry or clock')
    m=a[-1]; z=b[-1]; lam=a[0].ds.lam
    q0,q1=[evaluate(d,f,m,qcut) for d,f,_ in (a,b)]
    targets=[evaluate(d,f,pole,qcut)[3] for d,f,pole in (a,b)]
    base=1e-9*states[0]['x_end']
    if max(abs((pole*pole-base*base)/(lam*t/2)-1) for pole,t in zip((m,z),targets,strict=True))>3e-11:
        raise ArithmeticError('Independent gap equation failed')
    N0,S0,W0,J0=q0; N1,S1,W1,J1=q1
    # The cold shape can be isolated at the old cold count using the
    # normalized size of the regularized inverse-energy deficit.
    count=(N1-N0)*(1/m-S1/N1)
    shape=N0*(S0/N0-S1/N1)
    warm=W1-W0
    change=J1-J0
    if abs(count+shape+warm-change)>2e-15*max(J0,J1):
        raise ArithmeticError('Term-wise gap difference does not close')
    E=np.hypot(m,qcut)
    bound=N0*qcut*qcut/(m*E*(E+m))
    dm_bound=(lam/2)*abs(change)/(2*min(m,z))
    if abs(shape)>bound or abs(z-m)>dm_bound*(1+3e-7):
        raise ArithmeticError('Shape or monotone-root bound failed')
    return dict(evaluated_at_S2_mass=float(m),S4_mass=float(z),q_cold_upper=qcut,
                cold_counts=[N0,N1],cold_gap_fraction=(N0/m-S0)/J0,
                normalized_gap_rhs_differences=dict(cold_count=count/J0,cold_shape=shape/J0,
                                                    outside_cold=warm/J0,total=change/J0),
                worst_cold_shape_relative_rhs_bound=bound/J0,
                observed_fractional_mass_change=(z-m)/m,relative_mass_change_bound=dm_bound/m,
                scope='Algebraic comparison of the numerical nodal gap in floating point. No rate or other physical error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('states',nargs=2,type=Path)
    a=p.parse_args()
    print(json.dumps(run([read(f) for f in a.states]),indent=2,allow_nan=False))
