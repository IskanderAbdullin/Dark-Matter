"""Re-evaluate actual warm pressures at each chart's own numerical mass."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart, SPECS, load


def run(root,closure,mixed_lock,control):
    x,y,z=[p.read_bytes() for p in (closure,mixed_lock,control)]
    q,b=json.loads(x),json.loads(y)
    d=q['direct_comparison']
    if (q['existing_comparison_sha256'][0]!=hashlib.sha256(y).hexdigest()
            or d['budget_reference_sha256']!=hashlib.sha256(z).hexdigest()
            or len(d['contrasts'])!=5 or [d['charts'],b['charts']]!=[['mixed255','G8'],['mixed255','G7']]):
        raise ArithmeticError('Require the same three complete accepted histories')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    tag='2026-11-28'
    rows=[]
    for i,k in enumerate(('start','x01','x02','x04','x05')):
        m=d['contrasts'][i]
        x=m['x']
        common=m['fixed_native_mass']
        n=m['steps'][1]
        pins=[(b['uniform_trace'][i]['uncompressed_state_sha256'],b['contrasts'][i]['moments'][1]),
              (d['uniform_trace'][i]['uncompressed_state_sha256'],m['moments'][1]),
              (d['mixed_trace'][i]['state_sha256'],m['moments'][0])]
        ans=[]
        for g,(digest,old) in zip(('gauss7','gauss8','mixed255'),pins,strict=True):
            p=root/f'early-pole-fine-charge2-scalar4-{g}warm-{k}-state-{tag}.json.gz'
            card,h=read(p)
            if h!=digest:
                raise ArithmeticError('Underlying actual trajectory bytes changed')
            ctrl=control if g=='mixed255' else None
            a,_,_=gauss(card,base,common,x,expected_steps=n,budget_control=ctrl)
            if not np.array_equal(a,old):
                raise ArithmeticError('Common-mass moment did not reproduce exactly')
            _,chart,_=chart_for(card,budget_control=ctrl)
            native=chart.state(x,np.asarray(card['nodal_occupation']))[0]
            c,_,_=gauss(card,base,native,x,expected_steps=n,budget_control=ctrl)
            if a[0]!=c[0] or (c[1]-a[1])*(native-common)>0.:
                raise ArithmeticError('Mass remeasurement changed count or had the wrong pressure sign')
            ans.append(dict(chart=g,state_sha256=h,own_mass=native,pressure_common=a[1],pressure_native=c[1],
                            own_over_common_mass_minus_one=native/common-1.,pressure_native_over_common_minus_one=c[1]/a[1]-1.))
        v,w=[np.array([a[f'pressure_{s}'] for a in ans]) for s in ('common','native')]
        comp0=[v[1]/v[0]-1.,v[2]/v[1]-1.]
        comp1=[w[1]/w[0]-1.,w[2]/w[1]-1.]
        if np.max(abs(np.array(comp0)-comp1))>1e-9 or np.any(np.sign(comp0)!=np.sign(comp1)):
            raise ArithmeticError('Warm pressure ordering did not survive the measurement-mass change')
        rows.append(dict(x=x,common_measurement_mass=common,charts=ans,
                         contrasts_common_and_native=[comp0,comp1]))
    return dict(fields=['G8_over_G7_minus_one','mixed255_over_G8_minus_one'],rows=rows,
                scope='Remeasurement of existing warm-history states only, not re-evolution or cosmological pressure.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','closure','mixed-lock','control'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.closure,a.mixed_lock,a.control),indent=2,allow_nan=False))
