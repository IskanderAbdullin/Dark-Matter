"""Compare independently restarted charge-only nested trajectories on one scalar grid."""
import argparse
import json
from math import isfinite
from pathlib import Path


def load(p):
    return json.loads(p.read_text())['segments']


def summary(segments):
    if len(segments)!=4:
        raise ArithmeticError('A complete endpoint chain must have four restarts')
    grid=segments[0]['grid']
    div=segments[0]['subcells_per_archived_cell']
    counts=segments[0]['counts']
    rows=[]
    for s,start,end,steps in zip(segments,(.0075,.01,.02,.04),(.01,.02,.04,.05),(512,1024,1024,1024),strict=True):
        if (s['grid'],s['subcells_per_archived_cell'],s['counts'],s['starting_x'],s['x_end'],s['time_refinement'])!=(
                grid,div,counts,start,end,steps):
            raise ArithmeticError('Changed clock or chart in the chain')
        arms=s['checkpoint_trajectories']
        if len(arms)!=2 or any(len(r)!=8 for r in arms):
            raise ArithmeticError('Missing checkpoint samples on either time mesh')
        rows.append(dict(x=end,yield_hi=arms[1][-1]['scalar_yield'],
                         relative_time_drift=s['endpoint_yield_relative_drift']))
        if not (isfinite(rows[-1]['yield_hi']) and rows[-1]['yield_hi']>0 and
                0<=rows[-1]['relative_time_drift']<1e-3):
            raise ArithmeticError('Invalid yield or failed time comparison')
    return grid,div,counts,rows


def run(paths):
    tabs=[summary(load(p)) for p in paths]
    if len(tabs)<2:
        raise ValueError('Supply at least two complete nested continuations')
    g= tabs[0][0]
    for p,q in zip(tabs[:-1],tabs[1:],strict=True):
        if (p[0]!=g or q[0]!=g or q[1]!=2*p[1] or
            q[2]!=[p[2][0],2*p[2][1],2*p[2][2]]):
            raise ArithmeticError('Charge refinement must be nested at unchanged scalar order')
    lab=[z[1] for z in tabs]
    out=[]
    for rows in zip(*(z[3] for z in tabs),strict=True):
        ys=[d['yield_hi'] for d in rows]
        dif=[b-a for a,b in zip(ys[:-1],ys[1:],strict=True)]
        out.append(dict(x=rows[0]['x'],yields=ys,relative_time_drifts=[r['relative_time_drift'] for r in rows],
                        signed_relative_nested_changes=[b/a-1 for a,b in zip(ys[:-1],ys[1:],strict=True)],
                        signed_successive_increment_ratios=[b/a if a else None
                                                    for a,b in zip(dif[:-1],dif[1:],strict=True)]))
    return dict(grid=g,charge_divisions=lab,scalar_bins=tabs[0][2][0],matched_checkpoints=out,
                scope='Charge-only finite-grid comparison, not a limit or a physical transfer spectrum.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cards',nargs='+',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.cards),indent=2,allow_nan=False))
