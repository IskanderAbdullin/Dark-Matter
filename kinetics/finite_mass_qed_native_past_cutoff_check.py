"""Follow the accepted fine conditional source roots over exact past-band cutoffs."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from finite_mass_qed_native_scalar_bridge_join import field,normalize
from finite_mass_qed_native_scalar_bridge_receipt import read,run as receipt


def totals(low,rows,omega,mass):
    edges=np.r_[1e-6,.0003,[row['limits'][1] for row in rows]]
    chunks=[low]+[row['number_scalar_energy_and_centered_coefficients'] for row in rows]
    a=np.asarray([field(normalize(v,mass),omega,mass) for v in chunks],dtype=float)
    if min(np.ravel(a[:,:,:2]))<0 or np.any(a[:,:,2].sum(axis=0)<=0):
        raise ArithmeticError('Require nonnegative gross weights and positive whole slopes')
    c=np.array([[[fsum(a[:k,i,j]) for j in range(3)] for i in range(3)] for k in range(1,len(edges))])
    root=np.divide(c[:,:,1],c[:,:,2],out=np.full_like(c[:,:,1],np.nan),where=c[:,:,2]>0.)
    if not np.all(np.isfinite(root)) or np.any(root<=0.) or not np.all(np.diff(edges)>0.):
        raise ArithmeticError('A stored cutoff is not a valid positive conditional root')
    profile=np.r_[np.zeros((1,3,2)),c[:,:,:2]/c[-1,:,:2]]
    return edges,c,root,profile


def quantiles(edges,p,thresholds):
    if p[0]!=0 or p[-1]!=1 or np.any(np.diff(p)<0.):
        raise ArithmeticError('A positive gross prefix must span all normalized mass')
    k=np.searchsorted(p,thresholds)
    return np.c_[edges[k-1],edges[k]].tolist()


def run(root,low,cold,fixture):
    b=receipt(root,cold,'past-bridge',[0,1],'.gz',('fine',))
    if b!=fixture:
        raise ArithmeticError('Native fine bridge does not reproduce its independent ledger')
    b=b[0]
    omega,mass=b['omega'],b['m']
    a=next(x for x in json.loads(low.read_text()) if x['grid']=='fine')
    original=next(x for x in json.loads(cold.read_text()) if x['grid']=='fine')
    if (a['joined_support']!=[1e-6,.0003] or a['pieces'][0]['raw']!=original['integrated_number_scalar_energy_and_excess_coefficients']):
        raise ArithmeticError('Old low-state support or original frozen state changed')
    z=np.sum([normalize(r['raw'],mass) for r in a['pieces']],axis=0,dtype=np.longdouble)
    data=[totals(z,read(root/f'fine-past-bridge-{i}.jsonl.gz'),omega,mass) for i in (0,1)]
    x,y=data
    if not np.array_equal(x[0],y[0]):
        raise ArithmeticError('Paired supported prefix boundaries disagree')
    edges,c,r,p=y
    profiles_absolute_maximum=float(np.max(abs(x[3]-p)))
    relative_root_error=abs(x[2]/r-1)
    relative_root_max=float(np.max(relative_root_error))
    resolved=np.max(relative_root_error,axis=1)<=6e-5
    if profiles_absolute_maximum>3e-5 or not resolved[-1]:
        raise ArithmeticError(('Paired gross cutoff profiles or whole final root did not resolve',profiles_absolute_maximum,relative_root_max))
    Q=(.5,.9,.99)
    bounds=[[quantiles(edges,p[:,i,k],Q) for k in (0,1)] for i in range(3)]
    if bounds!=[[quantiles(edges,x[3][:,i,k],Q) for k in (0,1)] for i in range(3)]:
        raise ArithmeticError('A quantile bracket changed across meshes')
    out=[]
    for q in (.0005,.0008,.001,.0012,.0014):
        i=int(np.searchsorted(edges,q,'left'))
        out.append(dict(q=q,exact_bracketing_edges=edges[i-1:i+1].tolist(),gross_prefix_by_sector=p[i-1:i+1].tolist()))
    return dict(grid='fine',joined_support=[1e-6,.0015],channels=['scalar_number','scalar_energy','charged_energy'],
                separated_columns=['forward','inverse'],cutoff_boundaries=edges.tolist(),
                conditioned_roots_at_mesh_resolved_cutoffs=[v.tolist() if good else None for v,good in zip(r,resolved,strict=True)],
                maximum_relative_root_defect_by_cutoff=np.max(relative_root_error,axis=1).tolist(),
                unresolved_root_cutoff_count=int(np.count_nonzero(~resolved)),normalized_gross_prefix=p.tolist(),
                full_totals=c[-1].tolist(),quantile_levels=Q,quantile_boundaries_by_sector_and_gross_column=bounds,
                selected_cutoff_prefix_brackets=out,maximum_paired_absolute_cumulative_defect=profiles_absolute_maximum,
                maximum_paired_relative_root_defect=relative_root_max,
                scope='Prefix integrals of one frozen native past band; no within-panel interpolation or post-bridge inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','low','cold','fixture'):
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.low,a.cold,json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
