"""An external energy panel integrates the bare charged archive edge."""
from functools import lru_cache
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals
from finite_mass_qed_charged_closed_radial_pv_reference import full_radial
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes


def project(d,cut):
    M=d['x_end']
    m=M*d['mass_ratio']
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es,lo,hi,*_=intervals(M,s,m,Q)
    name,ma,shift=('parent',M,0.) if cut=='parent' else ('fermion',m,Es)
    edges=np.hypot(d[f'endpoint_{name}_momentum_edges'],ma)+shift
    occ=np.asarray(d[f'endpoint_{name}_occupations'])
    allowed=np.flatnonzero((edges[1:-1]>lo+4e-4)&(edges[1:-1]<hi-4e-4))
    k=1+allowed[np.argmax(abs(np.diff(occ)[allowed]))]
    center=edges[k]
    width=min(center-edges[k-1],edges[k+1]-center,center-lo,hi-center,0.024)/4
    p=np.sqrt((center-shift-ma)*(center-shift+ma))
    jump=occ[k]-occ[k-1]
    singular=solve_axes(M,m,s,center,cut)[1]*p*jump/2

    @lru_cache
    def evaluate(x):
        a=full_radial(d,center+x,cut,lo,hi,False)
        return np.array([a['step']['value'],a['linear_energy']['value']])

    def integrals(H,n,subtract):
        n,w=roots_legendre(n)
        g=(n+1)/2
        # Independently resolve the integrable outer log using either direct
        # dyadic panels or an analytic subtraction with a squared coordinate.
        if subtract:
            u=g*g
            W=.75*(1-u*u)*g*w
            pvals=np.array([evaluate(H*t)+evaluate(-H*t) for t in u])
            pvals[:,0]+=2*singular*np.log(H*u)
            ret=W@pvals
            ret[0]+=singular*(4/3-np.log(H))
        else:
            ret=np.zeros(2)
            stops=[0.,2.**-9,2.**-5,2.**-2,1.]
            for a,b in zip(stops[:-1],stops[1:],strict=True):
                u=a+(b-a)*g**2
                W=.75*(1-u*u)*(b-a)*g*w
                ret+=W@np.array([evaluate(H*t)+evaluate(-H*t) for t in u])
        return ret

    reports=[]
    for H in (width,width/4):
        a=integrals(H,6,True)
        b=integrals(H,10,True)
        c=integrals(H,6,False)
        if max(abs(a-b))>2e-3 or max(abs(c-b))>2e-3:
            raise ArithmeticError(f'Outer panel quadratures did not replay: {a,b,c}.')
        reports.append(dict(half_width=H,step_weighted_average=b[0],linear_weighted_average=b[1],
                            resolution_defect=float(max(abs(a-b))),direct_outer_replay_defect=float(max(abs(c-b)))))
    observed=reports[1]['step_weighted_average']-reports[0]['step_weighted_average']
    return dict(external_center_parent_energy=center,cut_species=cut,histogram_jump=float(jump),
                predicted_step_log_coefficient=singular,step_difference_under_fourfold_shrink=observed,
                leading_step_difference=singular*np.log(4.),panels=reports)


def run():
    results=[]
    for grid,cut in (('fine','daughter'),('third','parent')):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:
            data=json.load(f)
        results.append(dict(grid=grid,**project(data,cut)))
    return dict(scope='Bare charged cut with a normalized external parabolic window. No assignment of real-time cut factors or physical rate.',rows=results)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
