"""Check the late scalar-S4 clock at the two cold shape-factorial momenta."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_charge_cold_clock import optical
from early_pole_endpoint_charge_cold_clock_points import point,setup
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def run(coarse,fine,clock,reference):
    """One archived common .04 restart, using each clock's own final optical state."""
    a,sa=coarse
    b,sb=fine
    d,v=unpack(b,False,.05)
    keys=('grid','charge_div','scalar_div','cut','reference_x','counts','x_end')
    if a['steps']!=1024 or any(a[k]!=b[k] for k in keys) or a.get('warm_floor') is not None:
        raise ArithmeticError('Need the correctly labelled native-warm S4 clock')
    u=np.asarray(a['endpoint_occupation'])
    if u.shape!=(d.n,) or not np.all(np.isfinite(u)) or min(u)<0 or max(u[d.ns:])>1:
        raise ArithmeticError('Invalid recovered scalar state')
    if (clock['warm_div'],clock['x_start'],clock['x_end'],clock['steps'])!=(1,.04,.05,[1024,2048]):
        raise ArithmeticError('Expected the frozen S4 late trace comparison')
    c=clock['coarse_fine_count_pressure_and_warm_pressure']
    for i,y in enumerate((u,v)):
        if abs(fsum(d.ds.w*y[:d.ns])/c[i][0]-1)>3e-14:
            raise ArithmeticError('Archived states do not match the independent count receipt')
    if (reference['charge_subcells'],reference['grid'],reference['levels'],reference['state_sha256'][-1])!=(2,'fine',[1,2,4],sb):
        raise ArithmeticError('The independent pointwise factor swap uses a different state')
    params=[optical(d,y,.05) for y in (u,v)]
    opts=[setup(d,y,p) for y,p in zip((u,v),params,strict=True)]
    rows=[]
    for q in (3.5e-6,4.7e-6):
        r=next(s for s in reference['rows'] if s['q']==q and s['levels']==[2,4])
        t=next(s for s in r['rows'] if s['bits']==[1,1,1])
        A,B=t['bare_forward_inverse']
        pop=r['scalar_occupations'][-1]
        ref=np.array([(1+pop)*A,pop*B])
        (low,e),(high,f)=[point(d,w,o,q) for w,o in zip((u,v),opts,strict=True)]
        defect=float(np.max(abs(high[:2]/ref-1)))
        if defect>3e-10:
            raise ArithmeticError('Fine S4 control failed the independent factorial point receipt')
        rows.append(dict(q=q,coarse=low.tolist(),fine=high.tolist(),relative_coarse_fine=(low/high-1).tolist(),
                         independent_fine_point_max_relative_difference=defect,max_pair_mesh_defect=max(e,f)))
    return dict(grid='fine',x=.05,k=.0015,charge_subcells=2,cold_scalar_subcells=4,
                step_pair=[1024,2048],state_hashes_uncompressed=[sa,sb],optical_parameters=params,rows=rows,
                relative_count_coarse_fine=clock['coarse_fine_relative_change'][0],
                scope='Two conditional last-segment S4 sensitivities. Not full-band time or spatial convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('coarse','fine','clock','reference'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(read(a.coarse),read(a.fine),json.loads(a.clock.read_text()),json.loads(a.reference.read_text())),indent=2,allow_nan=False))
