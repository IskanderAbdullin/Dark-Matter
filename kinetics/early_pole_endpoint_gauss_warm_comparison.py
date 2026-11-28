"""Compare an evolved warm Gauss chart to the single-node warm trajectories."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_energy_moment_envelope import envelopes
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import fit,scaled_energy
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_spectral_moment_budget import integrate
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load
from early_pole_endpoint_warm_pressure_source_control import cell_weights
from early_pole_endpoint_warm_pressure_cells import partition


def chart_for(card, *, budget_control=None):
    schema=card['scalar_nodal_rule']
    if schema=='gauss-volume-v1' and 'warm_gauss_order' not in card:
        d=ScalarGauss2(load('fine')[0],SPECS['fine'],2,4,1e-5,gauss_floor=.01)
        return 2,d,d.scalar_repeats
    if schema=='gauss-volume-v2-general' and card.get('warm_gauss_order') in (3,4,5,6,7,8):
        from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
        order=card['warm_gauss_order']
        d=WarmGauss(load('fine')[0],order)
        return order,d,d.repeats
    if schema=='gauss-volume-v3-targeted':
        from early_pole_endpoint_targeted_warm_history import build
        d,meta=build()
        if any(card.get(k)!=v for k,v in meta.items()):
            raise ArithmeticError('Targeted warm archive does not bind the declared chart')
        return 'targeted-2xG4',d,d.repeats
    if schema=='gauss-volume-v4-fixed-budget' and budget_control is not None:
        from early_pole_endpoint_budgeted_warm_history import build
        d,meta=build(budget_control)
        if any(card.get(k)!=v for k,v in meta.items()):
            raise ArithmeticError('Fixed-budget warm archive does not bind the frozen selection')
        return 'mixed255',d,d.repeats
    raise ArithmeticError('Unknown or unbound warm Gauss node convention')


def gauss(card,base,m,x,*,expected_steps=None,budget_control=None):
    n=(1024 if x==.01 else 2048) if expected_steps is None else expected_steps
    if (card['grid'],card['x_end'],card['steps'],card['reference_x'],card['charge_div'],card['scalar_div'],
          card['cut'],card['gauss_floor'])!=('fine',x,n,X,2,4,1e-5,.01):
        raise ArithmeticError('Require the matched warm-only Gauss chart')
    order,d,repeats=chart_for(card,budget_control=budget_control)
    expected_cells=d.old_ns if order==2 else d.old_count
    if card['scalar_cell_count']!=expected_cells:
        raise ArithmeticError('Unrecognized native cell count')
    y=np.asarray(card['nodal_occupation'])
    if card['nodes']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or min(y)<0 or not np.all(np.isfinite(y)) or max(y[d.ns:])>1:
        raise ArithmeticError('Nodal archive mismatch')
    if not np.array_equal(d.ds.ed,base.ds.ed) or not np.array_equal(d.wd,base.wd) or not np.array_equal(d.wp,base.wp):
        raise ArithmeticError('Cold or charge geometry changed')
    P,W=[cell_weights(base.ds.ed,m,cut,64) for cut in (None,.01)]
    rows=[]
    cell_pairs=[]
    lower=[]
    upper=[]
    stable_width=[]
    max_fit_error=0.
    pos=np.r_[0,np.cumsum(repeats)]
    N=d.ds.w*y[:d.ns]
    for j,(a,b) in enumerate(zip(d.ds.ed[:-1],d.ds.ed[1:],strict=True)):
        ix=slice(pos[j],pos[j+1])
        u=fsum(N[ix])
        if not u:
            cell_pairs.append((0.,0.))
            continue
        if repeats[j]>1:
            s=fsum(N[ix]*scaled_energy(d.ds.p[ix],a,b,m))/u
            rec=fit(a,b,m,u,s)
            rho=u*rec['mean_pressure_numerator']
            e=np.array([a,b])
            exc=e*e/(np.hypot(e,m)+m)
            _,bounds,width=envelopes(exc[:1],exc[1:],exc[:1]+s*np.diff(exc),m)
            lo,hi=u*bounds[:,0,0]
            if not lo-4e-14*rho<=rho<=hi+4e-14*rho or a<.01:
                raise ArithmeticError('Warm fitted pressure escaped a conditional envelope')
            lower.append(lo)
            upper.append(hi)
            stable_width.append(u*m*m*width[0])
            max_fit_error=max(max_fit_error,rec['verify_error'])
            rows.append((u,rho,rho))
            cell_pairs.append((u,rho))
        else:
            rows.append((u,u*P[j],u*W[j]))
            cell_pairs.append((u,u*P[j]))
            lower.append(u*P[j])
            upper.append(u*P[j])
    z=np.array([fsum(col) for col in zip(*rows,strict=True)])
    if abs(z[0]/fsum(N)-1)>2e-14:
        raise ArithmeticError('Warm-only representative does not preserve total scalar count')
    lo,hi=fsum(lower),fsum(upper)
    w=fsum(stable_width)
    if abs((hi-lo)-w)>5e-16*z[1] or not lo-4e-14*z[1]<=z[1]<=hi+4e-14*z[1]:
        raise ArithmeticError('The conditionally allowable totals do not close')
    return z,dict(pressure_total_min_max=[lo,hi],relative_width_at_fixed_nodal_count_and_energy=w/z[1],
                  warm_exponential_fit_max_moment_error=max_fit_error),np.asarray(cell_pairs)


def run(s,ga,others,x=.05):
    d,y=unpack(s,False,x)
    m=d.state(x,y)[0]
    def direct(chart,y):
        a=integrate(chart,y,m,48)
        return np.array([fsum(a[:,0]),fsum(a[:,-1]),fsum(a[-2:,-1])])
    base=direct(d,y)
    G,envelope,cells=gauss(ga,d,m,x)
    weights=cell_weights(d.ds.ed,m,None,64)
    pop=d.ds.w*y[:d.ns]
    old=np.stack((pop,pop*weights),axis=-1)
    if abs(fsum(old[:,1])/base[1]-1)>2e-14 or abs(fsum(cells[:,1])/G[1]-1)>2e-14:
        raise ArithmeticError('Cell-attribution pressure check failed')
    delta=partition(old,cells)
    err=abs(fsum(delta[:,0])-(G[1]-base[1]))/base[1]
    if err>5e-15:
        raise ArithmeticError('Cell partition does not reproduce the accepted total')
    i=np.argsort(-abs(delta[:,0]))[:4]
    contribution=dict(components=('total_pressure_change','changed_cell_count','changed_within_cell_pressure_per_particle'),
                      total_relative_to_S4=[fsum(v)/base[1] for v in delta.T],
                      largest_absolute_changes=[dict(q=d.ds.ed[j:j+2].tolist(),relative_to_S4=(delta[j]/base[1]).tolist()) for j in i])
    out=[]
    for w in others:
        e,z=unpack(w,True,x)
        p=direct(e,z)
        out.append(dict(warm_div=w['warm_div'],moments=p.tolist(),relative_to_S4=(p/base-1).tolist(),
                        gauss_relative_to_this_warm=(G/p-1).tolist()))
    kind=ga['scalar_nodal_rule']
    tag=({} if kind=='gauss-volume-v1' else dict(warm_gauss_order=ga['warm_gauss_order'])
         if kind=='gauss-volume-v2-general' else dict(scalar_nodal_rule=kind,targeted_rule=ga['targeted_rule']))
    return dict(x=x,fixed_mass=float(m),**tag,columns=['count','pressure','pressure_q_ge_.01'],
                S4_reference=base.tolist(),warm_gauss_moments=G.tolist(),
                warm_gauss_conditional_pressure_envelope=envelope,warm_gauss_cell_attribution=contribution,
                warm_gauss_relative_to_S4=(G/base-1).tolist(),geometric_warm_charts=out,
                scope='Fixed-mass postprocessing; warm Gauss cells match nodal count and energy, other cells retain their constant convention. The interval conditions on evolved discrete moments, not solver error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('S4',type=Path)
    p.add_argument('gauss',type=Path)
    p.add_argument('geometric_warm',type=Path,nargs='*')
    p.add_argument('--x',choices=(.01,.02,.04,.05),default=.05,type=float)
    o=p.parse_args()
    print(json.dumps(run(read(o.S4),read(o.gauss),[read(p) for p in o.geometric_warm],o.x),indent=2,allow_nan=False))
