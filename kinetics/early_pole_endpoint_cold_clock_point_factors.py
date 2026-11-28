"""Attribute the final-clock effect at selected q with positive crossed controls."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_charge_cold_clock_points import setup,point
from early_pole_endpoint_cold_point_factorial import FACTORS,one
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,X,load


def run(states,clock):
    if clock['charts']!=[2,2,4,4] or clock['step_counts']!=[1024,2048,1024,2048]:
        raise ArithmeticError('The independent integral fixture or clocks changed')
    data=[]
    for n,i in ((2,0),(4,2)):
        d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,n,1e-5)
        recs=clock['full_state_panel_receipts'][i:i+2]
        g=[]; C=[]
        for j,((card,sha),ref) in enumerate(zip(states[i:i+2],recs,strict=True)):
            expected=('fine',.05,1024*(j+1),X,2,n,1e-5,[d.ns,d.nf,d.np])
            keys=('grid','x_end','steps','reference_x','charge_div','scalar_div','cut','counts')
            if tuple(card[k] for k in keys)!=expected or sha!=ref['state_sha256_uncompressed']:
                raise ArithmeticError('The matched band records do not bind these occupations')
            y=np.asarray(card['endpoint_occupation'])
            if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
                raise ArithmeticError('Invalid occupation')
            t=ref['diagnostic']['full_pole_and_snapshot']
            ch,cb,m,w,R=setup(d,y,[t[k] for k in ('scalar_mass','omega','residue')])
            g.append((m,(w,R),(ch,cb),d.ds.ed,y[:d.ns]))
            C.append((y,(ch,cb,m,w,R)))
        rows=[]
        for q in (3.5e-6,4.7e-6):
            v=one(g,q)
            A,B=[point(d,y,k,q)[0][:2] for y,k in C]
            decision=np.array(v.pop('fully_weighted_d2_over_d1'))
            check=float(np.max(abs((B/A)/decision-1)))
            if check>3e-8:
                raise ArithmeticError('Direct endpoint ratio does not reproduce the factorial ratio')
            paths=np.asarray(v['mean_log_factors_by_bare_chart'])
            scalar=np.asarray(v['scalar_population_log_factor'])
            log=np.log(decision)
            rows.append(dict(**v, fine_over_coarse=decision.tolist(), fine_over_coarse_log_share_by_factor=(np.vstack((paths,scalar))/log).tolist(),
                             direct_endpoint_ratio_defect=check))
        data.append(dict(scalar_subcells=n,coarse_and_fine_poles=[(g[j][0],*g[j][1]) for j in (0,1)],q_rows=rows))
    return dict(x=.05,k=.0015,columns=['gross_forward','gross_inverse'],
                fine_relative_to_coarse_factor_order=[*FACTORS,'scalar_occupation'],
                state_sha256_uncompressed=[p[1] for p in states],results=data,
                scope='Positive sampled factor replacements along both matched last-segment clocks. Neither crossed states nor a band-integrated separation of factors is implied.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--states',nargs=4,type=Path,required=True)
    a.add_argument('--matched',type=Path,required=True)
    p=a.parse_args()
    print(json.dumps(run([read(x) for x in p.states],json.loads(p.matched.read_text())),indent=2,allow_nan=False))
