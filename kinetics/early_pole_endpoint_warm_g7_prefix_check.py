"""Audit and attribute pressure differences between consecutive actual G6/G7 prefixes."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_chain_receipt import run as chain, EPOCHS, CLOCK
from early_pole_endpoint_gauss_warm_comparison import gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_pressure_cells import partition


def run(native, native_chain, reports6, states6, reports7, states7):
    n=len(native)
    if not 1<=n<=len(EPOCHS) or any(len(x)!=n for x in (reports6,states6,reports7,states7)):
        raise ArithmeticError('Need equal, consecutive, actual G6 and G7 prefixes')
    audits=[chain(p,q,prefix=True,order=k)['checked_epoch_pairs']
            for k,p,q in ((6,reports6,states6),(7,reports7,states7))]
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    out=[]
    for i,(x,c) in enumerate(zip(EPOCHS[:n],CLOCK[:n],strict=True)):
        b,bhash=read(native[i])
        q=native_chain['validated_chain'][i]
        y=np.asarray(b['endpoint_occupation'])
        if (bhash!=q['uncompressed_state_sha256'] or (b['x_end'],b['steps'])!=(x,2*c)
                or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1):
            raise ArithmeticError('Common accepted native chart does not match the epoch')
        m=d.state(x,y)[0]
        entries=[]
        for as_,s in zip(audits,(states6[i],states7[i]),strict=True):
            a,sh=read(s)
            if sh!=as_[i]['uncompressed_state_sha256']:
                raise ArithmeticError('A warm state changed after verifying the full prefix')
            entries.append(gauss(a,d,m,x,expected_steps=2*c))
        (a,e,f),(b,k,g)=entries
        lo0,hi0=e['pressure_total_min_max']
        lo1,hi1=k['pressure_total_min_max']
        t=b/a-1.
        interval=[lo1/hi0-1.,hi1/lo0-1.]
        z=partition(f,g)/a[1]
        whole=np.array([fsum(u) for u in z.T])
        hi=np.argsort(-abs(z[:,0]))[:5]
        if (abs(whole[0]-t[1])>2e-14 or abs(whole[0]-fsum(whole[1:]))>2e-14
                or not interval[0]-2e-14<=t[1]<=interval[1]+2e-14):
            raise ArithmeticError('Pressure change does not match its recomputed partition or envelope')
        out.append(dict(x=x,paired_steps=[c,2*c],common_native_sha256=bhash,mass=m,
                        G7_over_G6_minus_one=t.tolist(),pressure_interval=interval,
                        all_cell_partition=whole.tolist(),largest_native_cells=[
                            dict(interval_q=d.ds.ed[j:j+2].tolist(),contributions=z[j].tolist()) for j in hi]))
    return dict(fields=['count','pressure','warm_pressure'],
                partition_fields=['pressure_difference','changed_count','changed_pressure_per_particle'],
                rows=out,audit_receipts=audits,
                scope='Actual matched time pairs at a common native mass; count-energy pressure bounds are not a time-mesh or physical error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('native','reports6','states6','reports7','states7'):
        p.add_argument('--'+key,type=Path,nargs='+',required=True)
    p.add_argument('--native-chain',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,json.loads(a.native_chain.read_text()),a.reports6,a.states6,a.reports7,a.states7),
                     indent=2,allow_nan=False))
