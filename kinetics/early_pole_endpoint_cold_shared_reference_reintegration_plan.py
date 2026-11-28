"""Select one missing-support split per impossible coarser reference cell."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_panel_reweight_envelope import panel_excess
from early_pole_endpoint_cold_shared_reference_scan import feasibility
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels


def plan(left,right,scan):
    e=geometry(left)[0]
    r=geometry(right)[0]
    m=scan['fixed_reference_mass']
    ce=np.asarray(scan['coarse_parent_edges'])
    counts=np.asarray(scan['reference_count'])
    U=np.asarray(scan['reference_excess'])
    if (counts.shape!=U.shape or counts.shape!=(2,len(ce)-1) or np.any(counts<0) or np.any(U<0)):
        raise ArithmeticError('Invalid paired moment receipt')
    if np.max(abs(panel_excess(ce[:-1],ce[1:],m,48)*counts[0]-U[0])/np.maximum(U[0],1e-100))>8e-14:
        raise ArithmeticError('Coarser reference does not represent its original native constant state')
    mu=np.divide(U,counts,out=np.zeros_like(U),where=counts>0)
    failures=feasibility(e,ce,m,mu[1],counts[1])
    splits=[]
    for f in failures:
        a,b=f['parent_edges']
        i=int(np.searchsorted(e,a))
        j=int(np.searchsorted(e,b))
        if [e[i],e[j]]!=[a,b]:
            raise ArithmeticError('Original scalar boundary is absent')
        q=f['requested_mean']
        side=q>f['attainable_interval'][1]
        a,b=e[j-1:j+1] if side else e[i:i+2]
        for k in range(1,22):
            x=b-(b-a)*2**(-k) if side else a+(b-a)*2**(-k)
            if x in r:
                x=np.nextafter(x,(a if side else b))
            c,d=(x,b) if side else (a,x)
            mean=panel_excess(np.array([c]),np.array([d]),m,48)[0]
            if a<x<b and (mean>q if side else mean<q):
                splits.append([float(a),float(x),float(b)])
                break
        else:
            raise ArithmeticError('Single simple split did not recover allowed mean')
    q,changed=refine_panels(e,splits)
    if not np.array_equal(np.intersect1d(e,r),np.intersect1d(q,r)):
        raise ArithmeticError('A split unexpectedly changed the common display bins')
    for j in (0,1):
        if feasibility(q,ce,m,mu[j],counts[j]) or feasibility(r,ce,m,mu[j],counts[j]):
            raise ArithmeticError('Proposed split support does not fit both count-and-mean vectors')
    return dict(targeted_parent_panel_splits=splits,
                earlier_panel_count=len(e)-1,new_panel_count=len(q)-1,
                new_rates_to_integrate_at_each_mesh=len(changed),
                comparison_bin_count=len(np.intersect1d(q,r))-1,
                scope='Exact geometry and positive moment realizability only. None of the new rates is estimated or reused from a larger panel.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--left',type=Path,required=True)
    p.add_argument('--right',type=Path,required=True)
    p.add_argument('--scan',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(plan(a.left,a.right,json.loads(a.scan.read_text())),indent=2,allow_nan=False))
