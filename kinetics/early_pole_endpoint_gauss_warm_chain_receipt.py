"""Independently validate paired time logs and actual warm Gauss nodal states."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_chain_receipt import EPOCHS, read
from early_pole_endpoint_stencil_remap_pilot import SPECS, X, load
from early_pole_resolved_parent_reference import S

CLOCK = (2048, 512, 1024, 1024, 1024)


def verify(r, s, d, lo, hi, n, *, order=None):
    if order not in (None,3,4,5,6,7,8,9):
        raise ArithmeticError('Unsupported higher positive Gauss chart')
    cells_n, repeats = (d.old_ns,d.scalar_repeats) if order is None else (d.old_count,d.repeats)
    key = dict(grid='fine', charge_div=2, scalar_div=4, cut=1e-5, gauss_floor=.01,
               reference_x=X, nodes=[d.ns, d.nf, d.np], scalar_cell_count=cells_n,
               scalar_nodal_rule='gauss-volume-v1' if order is None else 'gauss-volume-v2-general')
    if order is not None:
        key['warm_gauss_order']=order
    return verify_chart(r,s,d,lo,hi,n,key,repeats)


def verify_chart(r,s,d,lo,hi,n,key,repeats):
    if any(r.get(k)!=v or s.get(k)!=v for k,v in key.items()):
        raise ArithmeticError('Warm Gauss chart is not the declared chart')
    if (r['starting_x'],r['x_end'],r['coarse_steps'],s['x_end'],s['steps'])!=(lo,hi,n,hi,2*n):
        raise ArithmeticError('Missing, reordered, or mislabelled clock interval')
    y = np.asarray(s['nodal_occupation'])
    if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid independent nodal state')
    tracks = r['checkpoint_trajectories']
    names = ('x','scalar_yield','pair_number_error','continuity_residual','max_charged','max_limited_extent')
    if len(tracks)!=2 or any(len(t)!=8 for t in tracks):
        raise ArithmeticError('Expected all paired cloud trace samples')
    u = np.asarray([[[row[k] for k in names] for row in t] for t in tracks])
    expected = lo*np.exp(np.log(hi/lo)*np.arange(1,9)/8)
    if (not np.all(np.isfinite(u)) or np.max(abs(u[:,:,0]/expected-1))>2e-14
            or np.any(u[:,:,1]<0) or np.min(u[:,:,2])<0 or np.max(u[:,:,2])>3e-13
            or np.max(abs(u[:,:,3]))>3e-10 or np.min(u[:,:,4])<0 or np.max(u[:,:,4])>1
            or np.min(u[:,:,5])<0 or np.max(u[:,:,5])>1):
        raise ArithmeticError('Trace clocks, continuity, bounds or paired numbers failed')
    v = u[1,-1,1]
    if v<=0:
        raise ArithmeticError('Terminal scalar yield not positive')
    drift = abs(u[0,-1,1]/v-1)
    w = (d.ds.w*y[:d.ns])/S
    agg = np.add.reduceat(w,np.r_[0,np.cumsum(repeats)[:-1]])
    cells = np.asarray(r['scalar_cell_yields'])
    mass = d.state(hi,y)[0]
    observed_drift=r['relative_time_drift']
    observed_mass=r['scalar_pole_mass']
    if (cells.shape!=agg.shape or not np.all(np.isfinite(cells)) or np.min(cells)<0
            or np.max(abs(cells-agg)/np.maximum(abs(agg),1e-105))>5e-14
            or abs(fsum(w)/v-1)>4e-14 or not np.isfinite(observed_drift) or abs(drift-observed_drift)>1e-14 or drift>1e-3
            or not np.isfinite(observed_mass) or observed_mass<=0 or abs(mass/observed_mass-1)>4e-14
            or abs(max(y[d.ns:])-u[1,-1,4])>4e-14):
        raise ArithmeticError('Saved state fails the reported counts, mass or time drift')
    return dict(x=hi,scalar_yield=v,pole_mass=mass,relative_time_drift=drift,
                maximum_absolute_continuity=float(np.max(abs(u[:,:,3]))),
                maximum_pair_number_error=float(np.max(u[:,:,2])))


def run(reports, states, *, prefix=False, order=None):
    k = len(reports)
    if k!=len(states) or not 1<=k<=len(EPOCHS) or not prefix and k!=len(EPOCHS):
        raise ArithmeticError('Need all epochs unless explicitly checking one chronological prefix')
    g=load('fine')[0]
    if order is None:
        d = ScalarGauss2(g,SPECS['fine'],2,4,1e-5,gauss_floor=.01)
    else:
        from early_pole_endpoint_warm_positive_gauss_ladder import WarmGauss
        d = WarmGauss(g,order)
    out = []
    previous=None
    for p,q,lo,hi,n in zip(reports,states,(X,*EPOCHS[:-1])[:k],EPOCHS[:k],CLOCK[:k],strict=True):
        r,rhash = read(p)
        s,shash = read(q)
        row = verify(r,s,d,lo,hi,n,order=order)
        if order is not None:
            if any('actual_previous_raw_sha256' not in c or c['actual_previous_raw_sha256']!=previous for c in (r,s)):
                raise ArithmeticError('Higher warm Gauss restart does not bind its previous fine bytes')
            row['verified_previous_raw_sha256']=previous
            previous=shash
        row.update(source_report_sha256=rhash,uncompressed_state_sha256=shash)
        out.append(row)
    x={} if order is None else dict(warm_gauss_order=order)
    return dict(complete=k==len(EPOCHS),scalar_nodal_rule='gauss-volume-v1' if order is None else 'gauss-volume-v2-general',nodal_counts=[d.ns,d.nf,d.np],**x,
                checked_epoch_pairs=out,
                scope='Archived trace and saved-state audit; hashes bind endpoint artifacts, not a claim of bytewise replay of every trajectory.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reports',nargs='+',type=Path,required=True)
    p.add_argument('--states',nargs='+',type=Path,required=True)
    p.add_argument('--prefix',action='store_true')
    p.add_argument('--order',type=int,choices=(3,4,5,6,7,8,9))
    a=p.parse_args()
    print(json.dumps(run(a.reports,a.states,prefix=a.prefix,order=a.order),indent=2,allow_nan=False))
