"""Test already selected warm node budgets on two additional frozen snapshots."""
import argparse
import hashlib
import json
from pathlib import Path

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_warm_budgeted_gauss_probe import evaluate


def run(selection,states,refs,comparisons):
    wire=selection.read_bytes()
    data=json.loads(wire)
    if ([u['x_end'] for u in states]!=[.02,.04] or [w['x'] for w in refs]!=[.02,.04]
        or (data['training_x'],data['held_check_x'])!=(.01,.05)
        or [u['maximum_total_scalar_nodes'] for u in data['choices']]!=[219,255]
        or any(w['warm_native_edges']!=data['warm_edges'] for w in refs)):
        raise ArithmeticError('Additional check must use the same selected rules on exact x=.02/.04 states')
    out=[]
    for card in data['choices']:
        path=card['warm_orders_by_cell']
        v=[evaluate(s,r,c,path) for s,r,c in zip(states,refs,comparisons,strict=True)]
        if any(row['nodes']!=card['applied_scalar_nodes'] for row in v):
            raise ArithmeticError('Held comparison changed the node budget')
        alt=[]
        for row in comparisons:
            alt.append([z['l1_signed_source_residual'] for z in row['charts'] if z['order'] in (6,7,8)])
        out.append(dict(warm_orders_by_cell=path,scalar_nodes=card['applied_scalar_nodes'],
                        check_results=v,uniform_G6_G7_G8_reference_L1_at_same_epochs=alt))
    return dict(training_x=.01,additional_held_x=[.02,.04],selection_sha256=hashlib.sha256(wire).hexdigest(),
                results=out,
                scope='Two further directly integrated frozen checks, not re-optimized, evolved or a statistical independence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True)
    for key in ('states','reference','comparisons'):
        p.add_argument('--'+key,type=Path,nargs=2,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.selection,*[[read(f) for f in getattr(a,k)] for k in ('states','reference','comparisons')]),indent=2,allow_nan=False))
