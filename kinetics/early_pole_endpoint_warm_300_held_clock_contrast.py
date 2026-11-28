"""Compare exact mixed300, G8 and G9 prefixes on the same S4 measurement clock."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import read, run as uniform
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_pressure_cells import partition
from early_pole_endpoint_warm_single_upgrade_history import audit


def run(z):
    k=len(z.native)
    truth=json.loads(z.native_chain.read_text())['validated_chain'][:k]
    if not 1<=k<=5 or any(len(s)!=k for s in (truth,z.g8_states,z.g8_reports,z.g9_states,z.g9_reports,z.mixed_states,z.mixed_reports)):
        raise ArithmeticError('Require matched actual and consecutive prefix lengths')
    x,y=[uniform(r,s,prefix=True,order=j)['checked_epoch_pairs'] for j,r,s in
            ((8,z.g8_reports,z.g8_states),(9,z.g9_reports,z.g9_states))]
    m=audit(z.control,z.mixed_reports,z.mixed_states)['prefix']
    a=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    out=[]
    for i,b in enumerate(truth):
        c,d=read(z.native[i])
        T,n=c['x_end'],c['steps']
        if d!=b['uncompressed_state_sha256'] or T!=b['x']:
            raise ArithmeticError('Underlying exact S4 reference changed')
        mass=a.state(T,np.asarray(c['endpoint_occupation']))[0]
        v=[];cells=[];e=[]
        pins=[x[i]['uncompressed_state_sha256'],y[i]['uncompressed_state_sha256'],m[i]['state_sha256']]
        for s,h,ctrl in zip((z.g8_states[i],z.g9_states[i],z.mixed_states[i]),pins,(None,None,z.control),strict=True):
            f,sh=read(s)
            if sh!=h:
                raise ArithmeticError('Remeasurement did not read the verified physical state bytes')
            value,disp,tab=gauss(f,a,mass,T,expected_steps=n,budget_control=ctrl)
            v.append(value);cells.append(tab);e.append(disp['pressure_total_min_max'])
        delta=np.asarray(v[1:])/v[0]-1.
        merges=[]
        for j in (1,2):
            q=partition(cells[0],cells[j])/v[0][1]
            p=[fsum(w) for w in q.T]
            lo,hi=e[j][0]/e[0][1]-1.,e[j][1]/e[0][0]-1.
            if abs(p[0]-delta[j-1,1])>2e-14 or not lo-2e-14<=p[0]<=hi+2e-14:
                raise ArithmeticError('Warm chart ledger disagreed with the independent count/pressure closure')
            merges.append(dict(all_cell_components=p,conditional_pressure_interval=[lo,hi]))
        out.append(dict(x=T,steps=n,common_native_mass=mass,states_sha256=pins,ratios_to_G8=delta.tolist(),details=merges))
    return dict(charts=['G8','G9','mixed300'],fields=['count','pressure','warm_pressure'],rows=out,
                scope='Verified actual prefix at common native measurement mass; no exactness, physical prediction, or adaptation during time evolution claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('control','native-chain'):
        p.add_argument('--'+name,type=Path,required=True)
    for name in ('native','g8-states','g8-reports','g9-states','g9-reports','mixed-states','mixed-reports'):
        p.add_argument('--'+name,type=Path,required=True,nargs='+')
    print(json.dumps(run(p.parse_args()),indent=2,allow_nan=False))
