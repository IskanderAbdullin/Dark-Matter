"""Select cold native midpoint evaluations by a guaranteed placement allowance."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import audit
from early_pole_endpoint_continuum_cold_band_pilot import refine_panels


def order(e,f,a,b):
    mass=[np.asarray(x,dtype=float)/fsum(x) for x in (a,b)]
    sizes=[x*np.diff(z) for x,z in zip(mass,(e,f),strict=True)]
    scores=sorted([(float(q),j,i) for j,v in enumerate(sizes) for i,q in enumerate(v)],reverse=True)
    return scores,fsum(q for q,_,_ in scores)


def run(fine,third,states,archive,fixture):
    a,b=[audit(p,s,archive,g) for p,s,g in zip((fine,third),states,('fine','third'),strict=True)]
    if [a[1],b[1]]!=fixture['states_sha256_uncompressed'] or [a[2],b[2]]!=fixture['gates']:
        raise ArithmeticError('Select from accepted matching native cold spectra only')
    U,V=[(a[0][i],b[0][i]) for i in (0,1)]
    e,f=V[0][0],V[1][0]
    if not np.array_equal(e,U[0][0]) or not np.array_equal(f,U[1][0]):
        raise ArithmeticError('Paired mesh topology changed')
    c=[order(s[0][0],s[1][0],s[0][1][:,0,0],s[1][1][:,0,0]) for s in (U,V)]
    (q,x),(r,y)=c
    if [v[1:] for v in q[:16]]!=[v[1:] for v in r[:16]] or abs(x-y)>2e-12:
        raise ArithmeticError('Selected bins are not stable between measured mesh pairs')
    rows=[]
    for n in (8,16):
        split=[]
        E=[]
        for j,z in enumerate((e,f)):
            centres=[[float(z[i]),float(z[i]+(z[i+1]-z[i])/2.),float(z[i+1])] for _,k,i in r[:n] if k==j]
            _,changed=refine_panels(z,centres) if centres else (z,set())
            split.append(centres)
            E.append(len(changed))
        gain=.5*fsum(g for g,_,_ in r[:n])
        rows.append(dict(original_panels_to_split=n,guaranteed_absolute_reduction_q=gain,
                         post_refinement_placement_allowance_q=y-gain,split_plans_by_grid=split,
                         new_coefficients_at_each_mesh_by_grid=E))
    return dict(grids=['fine','third'],states_sha256_uncompressed=fixture['states_sha256_uncompressed'],
                observable='forward_number',initial_placement_allowance_q=y,paired_total_allowance_difference_q=abs(x-y),
                plans=rows,scope='Targets only; no split rate or improved observed transport value is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('fine','third','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    for k in ('archive','fixture'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.fine,a.third,a.states,a.archive,json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
