"""Separate initialized hard-grid rate changes from finite first-window evolution."""
import argparse
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from early_pole_archived_shared_cold_window_reference import Cold


def measure(path,depth):
    with gzip.open(path,'rt') as f:A=json.load(f)['branches']
    s,t=A['steps'], A['cell-conservative']
    if any(v.get('aperture_subdivision',1)!=depth for v in (s,t)):
        raise ArithmeticError('Saved initializer has wrong hard geometry.')
    m0,m1=[np.asarray(p['snapshots'][0]['hard_moments']) for p in (s,t)]
    cold=np.asarray(s['snapshots'][0]['cold_bins'])
    ref=Cold(s['grid'],1,0,len(cold),True,exposure=s['duration'],
             reference_aperture_grid=s['aperture_grid'],archive_projection='steps',
             aperture_subdivision=depth)
    np.testing.assert_array_equal(ref.ns,cold)
    np.testing.assert_array_equal(ref.ns,t['snapshots'][0]['cold_bins'])
    np.testing.assert_array_equal(ref.M.init,m0)
    for a,b,c in zip(ref.M.edge_sets,s['hard_edge_sets'],t['hard_edge_sets'],strict=True):
        np.testing.assert_array_equal(a,b)
        np.testing.assert_array_equal(a,c)
    daughter=np.asarray(ref.M.edge_sets[2])
    low,high=1.4138888888888888,1.4666666666666666
    idx=np.flatnonzero((daughter[:-1]>=low-2e-14)&(daughter[1:]<=high+2e-14))
    if len(idx)!=2*depth or not np.all(np.diff(idx)==1):
        raise ArithmeticError('Nested aperture did not preserve the local charged interval.')
    np.testing.assert_allclose(daughter[[idx[0],idx[-1]+1]],[low,high],rtol=0.,atol=2e-14)
    i=ref.M.r[2].start+idx
    y=m0.copy()
    y[i]=m1[i]
    weights=np.r_[cold,0.,0.]
    n,p=[Cold.rhs(ref,np.r_[r.ravel(),weights],False)[ref.j0:-2] for r in (m0,y)]
    archived_rate=s['snapshots'][0].get('instantaneous_source_by_q')
    if archived_rate is not None:
        np.testing.assert_allclose(n,archived_rate,rtol=0.,atol=3e-14)
    δ=p-n
    return cold,dict(depth=depth,native_rate=fsum(n),subset_local_rate=fsum(p),
                     paired_initial_rate=fsum(δ),initial_pair_by_mode=δ.tolist(),
                     signed_modes=[int(np.count_nonzero(δ>0)),int(np.count_nonzero(δ<0))],
                     replaced_cells=[int(v) for v in i])


def run(paths,finitepath,local8):
    rows=[]
    with open(finitepath) as f:receipt=json.load(f)
    with gzip.open(local8,'rt') as f:raw8=json.load(f)
    q=np.array(raw8['cold_momentum_nodes_over_T'])
    for depth,file in zip((1,2,4,8),paths,strict=True):
        c,r=measure(file,depth)
        if rows:np.testing.assert_array_equal(anchor,c)
        else:anchor=c
        pt=next(s for s in receipt['levels'] if s['depth']==depth)
        gains=np.array(pt['terminal_mode_gains'])
        if len(gains)!=len(q) or np.any(gains<=0):
            raise ArithmeticError('Missing positive terminal increment.')
        tots=gains.reshape(4,8).sum(axis=1)
        r['terminal_gain_by_8_mode_band']=tots.tolist()
        r['terminal_fraction_in_highest_8_modes']=float(tots[-1]/tots.sum())
        r['paired_terminal_rms_momentum_over_T']=float(np.sqrt(np.dot(q*q,gains)/gains.sum()))
        r['paired_first_window_mean']=pt['paired_mean']
        r['negative_to_nonnegative_brackets']=pt['first_cross_brackets']
        r['lowest_paired_count']=min(pt['paired_counts'])
        r['exposure_of_lowest_paired_count']=receipt['exposure'][int(np.argmin(pt['paired_counts']))]
        rows.append(r)
    latest=rows[-1]
    match=abs(latest['subset_local_rate']-raw8['local_archive_support']['initialized_cold_rate'])
    if match>2e-15:
        raise ArithmeticError('Initialized rate disagrees with the actual depth-eight run receipt.')
    return dict(scope='Matched cold32 static initialization and observed finite first-window responses only.',
                rows=rows,initial_rate_depth_adjacent_changes=np.diff([r['paired_initial_rate'] for r in rows]).tolist(),
                integrated_mean_depth_adjacent_changes=np.diff([r['paired_first_window_mean'] for r in rows]).tolist(),
                depth8_actual_initialized_rate_difference=match)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('initializers',nargs=4)
    p.add_argument('--finite',required=True)
    p.add_argument('--local8',required=True)
    a=p.parse_args()
    print(json.dumps(run(a.initializers,a.finite,a.local8),indent=2,allow_nan=False))
