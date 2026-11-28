"""Direct mixed255/G8 continuation with transitive checks against both G7 comparisons."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_warm_fixed_budget_contrast import run as compare


def verify(d,p,q):
    if (d['charts']!=['mixed255','G8'] or p['charts']!=['mixed255','G7']
            or any(len(x)!=5 for x in (d['contrasts'],p['contrasts'],q['rows']))
            or d['budget_reference_sha256']!=p['budget_reference_sha256']
            or d['mixed_trace']!=p['mixed_trace']
            or [p['uniform_trace'],d['uniform_trace']]!=q['audit_receipts']):
        raise ArithmeticError('Require all three exact complete histories')
    out=[]
    for a,b,c in zip(d['contrasts'],p['contrasts'],q['rows'],strict=True):
        if ((a['x'],a['steps'],a['common_native_sha256'],a['fixed_native_mass'])
                !=(b['x'],b['steps'],b['common_native_sha256'],b['fixed_native_mass'])
                or (a['x'],a['steps'],a['common_native_sha256'],a['fixed_native_mass'])
                !=(c['x'],c['paired_steps'],c['common_native_sha256'],c['mass'])
                or a['moments'][0]!=b['moments'][0]):
            raise ArithmeticError('All measured moments must share the exact clock, mass and mixed bytes')
        actual=np.asarray(a['mixed255_over_G8_minus_one'])
        z=(1.+np.asarray(b['mixed255_over_G7_minus_one']))/(1.+np.asarray(c['G8_over_G7_minus_one']))-1.
        err=np.max(abs(actual-z))
        if err>2e-14:
            raise ArithmeticError('Direct three-chart measurements do not close')
        lo,hi=a['conditional_pressure_contrast']
        out.append(dict(x=a['x'],measurement_mass=a['fixed_native_mass'],maximum_triangle_defect=err,
                        direct_mixed255_over_G8_minus_one=actual.tolist(),conditional_pressure_interval=[lo,hi],
                        conditional_pressure_order='mixed_lower' if hi<0. else 'mixed_higher' if lo>0. else 'unresolved',
                        pressure_cell_partition=a['cell_partition']))
    return out


def run(args):
    p,q=[x.read_bytes() for x in (args.mixed_g7_lock,args.g8_g7_lock)]
    z=compare(args.control,json.loads(args.native_chain.read_text()),args.native_states,args.mixed_reports,args.mixed_states,
              args.g8_reports,args.g8_states,attribution=True,uniform_order=8)
    out=verify(z,json.loads(p),json.loads(q))
    return dict(direct_comparison=z,verified_three_chart_closure=out,existing_comparison_sha256=[hashlib.sha256(x).hexdigest() for x in (p,q)],
                scope='Actual independently evolved fixed spatial layouts at matched clocks and mass. Conditional count-energy envelope is not time error, a continuum limit, or accuracy ranking.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('control','native-chain','mixed-g7-lock','g8-g7-lock'):
        p.add_argument('--'+k,type=Path,required=True)
    for k in ('native-states','mixed-reports','mixed-states','g8-reports','g8-states'):
        p.add_argument('--'+k,type=Path,nargs='+',required=True)
    print(json.dumps(run(p.parse_args()),indent=2,allow_nan=False))
