"""Compare recovered original coarse pressure states with the accepted G6/G7 continuation."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import verify
from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_resolved_parent_reference import S


def run(native, prefix, raw_reports, starts, coarse, fine, coarse_replays, *, terminal=False):
    lo,x,n=(.04,.05,1024) if terminal else (.0075,.01,512)
    card,h=read(native)
    d,y=unpack(card,False,x)
    m=d.state(x,y)[0]
    r=prefix['rows'][-1]
    if (r['x'],r['common_native_sha256'],r['mass'],r['paired_steps'])!=(x,h,m,[n,2*n]):
        raise ArithmeticError('Pressure clock is not fixed at the independently paired native mass')
    blocks=[]
    for i,k in enumerate((6,7)):
        report,rh=read(raw_reports[i])
        start,sh=read(starts[i])
        c,ch=read(coarse[i])
        f,fh=read(fine[i])
        p,ph=read(coarse_replays[i])
        left,right=prefix['audit_receipts'][i][-2:]
        if (sh,rh,fh)!=(left['uncompressed_state_sha256'],right['source_report_sha256'],right['uncompressed_state_sha256']):
            raise ArithmeticError('Replayed states do not bind the accepted prefix')
        if p['x_start']!=start['x_end'] or (p['x_start'],p['x_end'],p['steps'],p['warm_gauss_order'],
               p['fixed_pressure_mass'])!=(lo,x,[n,2*n],k,m) or not 0<=p['maximum_trace_absolute_difference']<=5e-13:
            raise ArithmeticError('Coarse replay trace failed time or mass acceptance')
        order,chart,_=chart_for(f)
        if order!=k:
            raise ArithmeticError('Pressure order mismatch')
        verify(report,f,chart,lo,x,n,order=k)
        if any(a['actual_previous_raw_sha256']!=sh for a in (c,f,report)) or chart_for(c)[0]!=k:
            raise ArithmeticError('Original restart hash or coarse order changed')
        Q=[gauss(a,d,m,x,expected_steps=z)[0] for a,z in ((c,n),(f,2*n))]
        trace=report['checkpoint_trajectories'][0][-1]['scalar_yield']
        dy=fsum(chart.ds.w*np.asarray(c['nodal_occupation'][:chart.ns]))/S
        if abs(dy/trace-1)>4e-14 or np.max(abs(Q/np.array(
                [p['coarse_count_and_pressures'],p['fine_count_and_pressures']])-1.))>5e-15:
            raise ArithmeticError('Coarse state failed independent yield or pressure replay')
        blocks.append(dict(order=k,source_report_sha256=rh,coarse_replay_sha256=ph,
                           coarse_state_sha256=ch,fine_state_sha256=fh,
                           maximum_trace_abs_difference=p['maximum_trace_absolute_difference'],
                           values=[v.tolist() for v in Q]))
    v=np.array([q['values'] for q in blocks])
    c=v[1]/v[0]-1.
    if np.max(abs(c[1]-r['G7_over_G6_minus_one']))>5e-15:
        raise ArithmeticError('Pressure clock table does not reproduce the independent fine contrast')
    return dict(x=x,starts=lo,steps=[n,2*n],measurement_mass=m,
                fields=['count','pressure','warm_pressure'],charts=blocks,
                G7_over_G6_minus_one=c.tolist(),fine_minus_coarse_contrast=(c[1]-c[0]).tolist(),
                scope=('Conditional replay of the terminal continuation at each order, retaining its own x=.04 state. No full-history clock bound.'
                       if terminal else 'Conditional replay of the first continuation at each order, retaining its own startup state. No startup or full-history clock bound.'))


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--native',type=Path,required=True)
    a.add_argument('--prefix',type=Path,required=True)
    for name in ('raw-reports','starts','coarse','fine','coarse-replays'):
        a.add_argument('--'+name,type=Path,nargs=2,required=True)
    a.add_argument('--terminal',action='store_true')
    p=a.parse_args()
    print(json.dumps(run(p.native,read(p.prefix)[0],p.raw_reports,p.starts,p.coarse,p.fine,p.coarse_replays,terminal=p.terminal),
                     indent=2,allow_nan=False))
