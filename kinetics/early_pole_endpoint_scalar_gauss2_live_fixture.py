"""Test a genuine two-node scalar evolution for short fixed-support controls only."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import expit

from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_scalar_positive_cell_quadrature import volume_rule
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_resolved_parent_reference import S


class ScalarGauss2(ScalarSubshellChart):
    def __init__(self,a,shape,charge_div,scalar_div,cut,*,gauss_floor=None):
        super().__init__(a,shape,charge_div,scalar_div,cut)
        q,w,_=volume_rule(self.ds.ed[:-1],self.ds.ed[1:])
        self.old_ns=self.ns
        correction=self.ds.w/w.sum(axis=1)
        self.volume_roundoff=float(np.max(abs(correction-1)))
        if self.volume_roundoff>1e-10:
            raise ArithmeticError('The original count uses a materially different shell volume')
        # Keep the stored scalar cell volume (legacy difference of cubes);
        # assign its two parts by the positive analytic normalized rule.
        w*=correction[:,None]
        if not np.allclose(w.sum(axis=1),self.ds.w,rtol=3e-15,atol=0):
            raise ArithmeticError('Gauss2 changed an incoming scalar cell count')
        if gauss_floor is None:
            take=np.ones(self.old_ns,dtype=bool)
        else:
            if gauss_floor<cut or not np.isfinite(gauss_floor):
                raise ValueError('An optional Gauss floor must leave the cold shells unchanged')
            take=self.ds.ed[:-1]>=gauss_floor
        if not np.any(take):
            raise ValueError('No selected scalar shell for Gauss nodes')
        sizes=np.where(take,2,1)
        self.scalar_repeats=sizes
        self.gauss_shells=int(np.count_nonzero(take))
        select=np.arange(2)[None,:]<sizes[:,None]
        q[~take,0]=self.ds.p[~take]
        w[~take,0]=self.ds.w[~take]
        self.ds.p=q[select]
        self.ds.w=w[select]
        self.ds.n=len(self.ds.w)
        self.ns=self.ds.n
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]
        self.initial=self.transfer_step(self.initial)

    def transfer_step(self,state):
        if state.shape!=(self.old_ns+self.nf+self.np,):
            raise ArithmeticError('Expected one constant occupation per original scalar cell')
        return np.r_[np.repeat(state[:self.old_ns],self.scalar_repeats),state[self.old_ns:]]


def probe(d,f,x):
    m,E,_,_=d.state(x,f)
    raw,projection=d.chart(x,m,8)
    thermal=np.r_[1/np.expm1(E[:d.ns]),expit(-E[d.ns:])]
    _,F,B=event(raw,thermal,d)
    eq=float(np.max(abs(F-B)/np.maximum(F,1e-99)))
    event_rates,_,_=event(raw,f,d)
    number=[fsum(t) for t in np.split(event_rates,(d.ns,d.iP))]
    pair_err=max(abs(number[0]-number[1]),abs(number[0]+number[2]))
    energy=abs(fsum(E*event_rates))
    if eq>3e-13 or pair_err>2e-13 or energy>2e-13:
        raise ArithmeticError('Live Gauss2 chart failed moment or equilibrium test')
    dt=5e-5
    end=x*np.exp(dt)
    arms=[drive(d,end,n,x_start=x,initial=f) for n in (32,64)]
    initial=fsum(d.w[:d.ns]*f[:d.ns])/S
    out=[]
    for n,(trace,y) in zip((32,64),arms,strict=True):
        if len(trace)!=8 or min(y)<0 or max(y[d.ns:])>1 or any(
                t['pair_number_error']>3e-13 or abs(t['continuity_residual'])>3e-13 for t in trace):
            raise ArithmeticError('Short live chart failed explicit bounded continuation')
        out.append(dict(steps=n,scalar_number_gain=trace[-1]['scalar_yield']-initial,
                        max_pair_error=max(t['pair_number_error'] for t in trace),
                        max_continuity=max(abs(t['continuity_residual']) for t in trace)))
    if abs(out[1]['scalar_number_gain']/out[0]['scalar_number_gain']-1)>3e-4:
        raise ArithmeticError('Time sensitivity masks this short live chart contrast')
    return dict(scalar_nodes=d.ns,initial_gap=float(m),thermal_affine_relative_defect=eq,
                event_max_energy_projection=projection,event_scalar_gain=number[0],
                event_pair_absolute_defect=pair_err,event_energy_absolute_defect=energy,
                short_live_arms=out)


def run(card,gauss_floor=None):
    x=card['x_end']; g=card['grid']
    t=dict(charge_div=card['charge_div'],scalar_div=card['scalar_div'],cut=card['cut'])
    if x!=.05 or g!='fine' or card['reference_x']!=X or card.get('warm_floor') is not None:
        raise ArithmeticError('Use an unmodified warm-grid frozen fine endpoint')
    a=load(g)[0]
    old=ScalarSubshellChart(a,SPECS[g],**t)
    high=ScalarGauss2(a,SPECS[g],**t,gauss_floor=gauss_floor)
    y=np.asarray(card['endpoint_occupation'])
    q=high.transfer_step(y)
    if card['counts']!=[old.ns,old.nf,old.np] or not np.all(np.isfinite(y)) or min(y)<0 or max(y[old.ns:])>1:
        raise ArithmeticError('State does not match declared input chart')
    if abs(fsum(old.w[:old.ns]*y[:old.ns])/fsum(high.w[:high.ns]*q[:high.ns])-1)>3e-14:
        raise ArithmeticError('Gauss2 lost number from same reconstructed initial shape')
    if any(not np.array_equal(getattr(old,k),getattr(high,k)) for k in ('qd','qp','wd','wp')):
        raise ArithmeticError('Gauss2 altered the charged chart')
    p,r=[probe(d,f,x) for d,f in ((old,y),(high,q))]
    return dict(grid=g,x=x,chart=t,gauss_floor=gauss_floor,gauss_shells=high.gauss_shells,old_single=p,new_gauss2=r,
                max_legacy_volume_roundoff_renormalization=high.volume_roundoff,
                relative_change_in_fine_short_gain=r['short_live_arms'][-1]['scalar_number_gain']/
                    p['short_live_arms'][-1]['scalar_number_gain']-1,
                scope='Independent scalar nodal occupations in a short local history; no whole trajectory or physical distribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('--gauss-floor',type=float)
    a=p.parse_args()
    print(json.dumps(run(read(a.state),a.gauss_floor),indent=2,allow_nan=False))
