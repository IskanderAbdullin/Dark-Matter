"""Validate a chain of identical endpoint-chart restarts."""
import argparse
import json
from pathlib import Path


def run(paths):
    parts=[]
    for p in paths:
        x=json.loads(Path(p).read_text())
        report=x['checkpoint_trajectories']
        if len(report)!=2 or x['endpoint_yield_relative_drift']>1e-3:
            raise ArithmeticError('Both declared temporal meshes are required')
        if any(d['pair_number_error']>3e-13 or d['max_charged']>1
               or abs(d['continuity_residual'])>3e-10 for r in report for d in r):
            raise ArithmeticError('A reported sample did not pass the declared budgets')
        tail=report[-1][-1]
        if abs(tail['x']/x['x_end']-1)>5e-15:
            raise ArithmeticError('Stage ended at a different x')
        if parts:
            prev=parts[-1]
            if x['grid']!=prev['grid'] or x['counts']!=prev['counts'] or x['subcells_per_archived_cell']!=prev['subcells_per_archived_cell']:
                raise ArithmeticError('Changed architecture during continuation')
            if x['starting_x']!=prev['x_end'] or abs(x['starting_scalar_yield']/prev['checkpoint_trajectories'][-1][-1]['scalar_yield']-1)>3e-15:
                raise ArithmeticError('Restart yield or clock mismatch')
        if x.get('original_archive_yield'):
            x['signed_endpoint_yield_change_from_original']=tail['scalar_yield']/x['original_archive_yield']-1
        parts.append(x)
    return dict(scope='Matched scalar/number and energy checks on one fixed spatial chart, not spatial convergence',
                segments=parts)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('segment_reports',nargs='+')
    print(json.dumps(run(p.parse_args().segment_reports),indent=2,allow_nan=False))
