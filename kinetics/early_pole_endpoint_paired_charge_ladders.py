"""Pair entire native charge-resolution ladders without combining scalar grids."""
import argparse
import json
from math import isfinite
from pathlib import Path


def run(a,b):
    if a['grid']!='fine' or b['grid']!='third' or a['charge_divisions']!=b['charge_divisions']:
        raise ArithmeticError('Need fine and third at identical charge increments')
    if a['scalar_bins']==b['scalar_bins'] or len(a['matched_checkpoints'])!=4 or len(b['matched_checkpoints'])!=4:
        raise ArithmeticError('Need both complete independent scalar meshes')
    out=[]
    for u,v,x in zip(a['matched_checkpoints'],b['matched_checkpoints'],(.01,.02,.04,.05),strict=True):
        n=a['charge_divisions']
        if u['x']!=x or v['x']!=x or len(u['yields'])!=len(n) or len(v['yields'])!=len(n):
            raise ArithmeticError('Mismatched endpoint times or subdivisions')
        f,t=u['yields'],v['yields']
        if any(not isfinite(y) or y<=0 for y in f+t):
            raise ArithmeticError('Invalid scalar number estimate')
        rel=[z/y-1 for y,z in zip(f,t,strict=True)]
        out.append(dict(x=x,fine=f,third=t,third_relative_to_fine=rel,
                        fine_nested_changes=u['signed_relative_nested_changes'],
                        third_nested_changes=v['signed_relative_nested_changes']))
    return dict(charge_subcells=n,scalar_counts=[a['scalar_bins'],b['scalar_bins']],
                matched_x_comparisons=out,scope='End-to-end finite chart comparison; smaller peer spread is not convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('fine',type=Path)
    p.add_argument('third',type=Path)
    a=p.parse_args()
    print(json.dumps(run(json.loads(a.fine.read_text()),json.loads(a.third.read_text())),indent=2,allow_nan=False))
