"""Measure propagated prior-clock changes at two cold points without re-timing the last interval."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_charge_cold_clock_points import point,setup
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(paths,reports):
    archive=[read(p) for p in paths]
    out=[]
    for n in (2,4):
        i=3*((n-2)//2)
        (old,sh),(a,ah),(b,bh)=archive[i:i+3]
        z=reports[(n-2)//2]
        d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
        k=dict(grid='fine',charge_div=2,scalar_div=n,cut=1e-5,counts=[d.ns,d.nf,d.np],reference_x=X)
        if any(any(x.get(j)!=v for j,v in k.items()) for x in (old,a,b)):
            raise ArithmeticError('Unexpected fixed scalar topology')
        if (old['x_end'],old['steps'],a['x_end'],a['steps'],b['x_end'],b['steps'])!=(.04,1024,.05,2048,.05,2048):
            raise ArithmeticError('Earlier clock changed the final clock or support')
        if (a['actual_prior_state_sha256']!=sh or z['inputs_sha256'][2]!=bh or
            z['alternate_path']!=[1024,2048] or z['original_path']!=[2048,2048] or z['archive_prior_coarse_trace_defect']>5e-13):
            raise ArithmeticError('Does not bind actual hashed earlier restart and true archived fine terminal')
        v=[np.asarray(x['endpoint_occupation']) for x in (a,b)]
        if any(x.shape!=(d.n,) or not np.all(np.isfinite(x)) or min(x)<0 or max(x[d.ns:])>1 for x in v):
            raise ArithmeticError('Invalid archived positivity')
        e=d.ds.ed
        lo,hi=np.maximum(e[:-1],1e-6),np.minimum(e[1:],1.5e-5)
        u=np.maximum(0.,hi-lo)*(lo*lo+lo*hi+hi*hi)/3
        f,g=[u*x[:d.ns] for x in v]
        tv=fsum(abs(f/fsum(f)-g/fsum(g)))/2
        o=[setup(d,y,p) for y,p in zip(v,z['final_alternate_and_original_poles'],strict=True)]
        result=[]
        for q in (3.5e-6,4.7e-6):
            (R,c),(S,t)=[point(d,y,h,q) for y,h in zip(v,o,strict=True)]
            result.append(dict(q=q,alternate_over_original_minus_one=(R/S-1).tolist(),maximum_paired_mesh_defect=max(c,t)))
        out.append(dict(scalar_subcells=n,original_and_alternate_state_sha256=[bh,ah],
                        original_coarse_recovery_state_sha256=sh,scalar_band_normalized_TV=tv,
                        scalar_band_count_alternate_over_original_minus_one=fsum(f)/fsum(g)-1,rows=result))
    return dict(x=.05,k=.0015,conditional_coefficients=['forward_number','inverse_number','affine_slope'],
                results=out,scope='Same terminal clock with prior .02–.04 timing changed. Point-only response, not a q-integrated or complete history error estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--states',nargs=6,type=Path,required=True,metavar=('S2_X04','S2_ALT','S2_REF','S4_X04','S4_ALT','S4_REF'))
    p.add_argument('--reports',nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.states,[json.loads(p.read_text()) for p in a.reports]),indent=2,allow_nan=False))
