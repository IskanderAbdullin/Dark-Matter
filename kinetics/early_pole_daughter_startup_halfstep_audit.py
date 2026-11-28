"""Instrument only the first shared outer interval without changing solver inputs."""
import argparse
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

import early_pole_cold_screened_split_reference as split
from early_pole_local_daughter_elastic_subcycle_reference import run as evolve


def capture(bundle, case, n, k):
    with gzip.open(bundle,'rt') as f:
        modes=len(json.load(f)['branches']['steps']['snapshots'][0]['cold_bins'])
    if modes!=8 or n*k!=64:
        raise ValueError('Use the original fixed-inner eight-mode experiment.')
    ivp=split.solve_ivp
    runner=split.run
    i=n//16
    steps=[]
    with TemporaryDirectory() as d:
        path=Path(d)/'partial.json'
        def bounded(*args,**kwargs):
            budget=kwargs['elastic_energy_budget']
            v=runner(*args,**kwargs,checkpoint=path,stop_after=i,checkpoint_every=i)
            return dict(v,elastic_energy_budget=budget)
        def audit(fun,span,seed,**opts):
            x=np.asarray(seed)[-modes-2:-2].copy()
            y=ivp(fun,span,seed,**opts)
            z=y.y[-modes-2:-2,-1].copy()
            if not y.success:
                raise ValueError('Invalid slow solve.')
            steps.append(dict(width=span[1]-span[0],start=x.tolist(),end=z.tolist(),
                              mean_rate=float((z-x).sum()/(span[1]-span[0]))))
            return y
        with patch.object(split,'solve_ivp',audit),patch.object(split,'run',bounded):
            result=evolve(bundle,case,n,k)
        rec=json.loads(path.read_text())
    if not result['checkpoint_only'] or rec['steps_done']!=i or len(steps)!=2*i:
        raise ValueError('Wrong startup window.')
    for a,b in zip(steps[:-1],steps[1:],strict=True):
        np.testing.assert_array_equal(a['end'],b['start'])
    np.testing.assert_allclose(sum(steps[-1]['end'])-sum(steps[0]['start']),
                               result['scalar_gain_so_far'],rtol=0.,atol=2e-16)
    peak=max(abs(q['rho_drift']) for v in rec['controls'] for q in v['elastic_substeps'])
    if peak>=result['elastic_substep_energy_budget']:
        raise ValueError('Original per-inner-step budget failed.')
    return dict(case=case,macro=n,inner=k,half_steps=steps,
                gain=result['scalar_gain_so_far'],peak_abs_inner_energy_defect=peak)


def compare(bundle, originals):
    out={}
    for key,n,k in [('outer16',16,4),('outer64',64,1)]:
        a,b=[capture(bundle,x,n,k) for x in ('native','local-daughter')]
        for x,v in zip((a,b),originals[key],strict=True):
            idx=n//16
            expected=v['state_history'][idx]-v['state_history'][0]
            np.testing.assert_array_equal(x['gain'],expected)
        d=[]
        for x,y in zip(a['half_steps'],b['half_steps'],strict=True):
            nx=np.asarray(x['end'])-x['start']
            ny=np.asarray(y['end'])-y['start']
            d.append(float(np.sum(ny-nx)))
        out[key]=dict(native=a,local=b,halfstep_paired_count_gains=d,
                      total_paired_count_gain=b['gain']-a['gain'])
    return dict(cases=out,outer16_minus_outer64=float(
                   out['outer16']['total_paired_count_gain']-out['outer64']['total_paired_count_gain']),
                scope='Instrumented first interval only; prior checkpoint reproduces exactly.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    for k in ('outer16','outer64'):
        p.add_argument('--'+k,nargs=2,required=True,help='Original native and local states.')
    a=p.parse_args()
    print(json.dumps(compare(a.bundle,{k:[json.loads(Path(f).read_text()) for f in v]
                     for k,v in vars(a).items() if k!='bundle'}),indent=2,allow_nan=False))
