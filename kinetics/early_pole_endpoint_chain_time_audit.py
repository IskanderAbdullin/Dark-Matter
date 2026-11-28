"""Compare time ladders for the whole fine endpoint-cell trajectory."""
import argparse
import json
from pathlib import Path


def run(old,high):
    base=json.loads(old.read_text())['segments']
    ref=[json.loads(f.read_text()) for f in high]
    if len(base)!=len(ref) or [z['x_end'] for z in ref] != [.01,.02,.04,.05]:
        raise ArithmeticError('Complete equal-epoch trajectories required')
    rows=[]
    prev=None
    for a,b in zip(base,ref,strict=True):
        if any(a[k]!=b[k] for k in ('grid','subcells_per_archived_cell','counts','starting_x','x_end')):
            raise ArithmeticError('Changed spatial chart or epoch during time refinement')
        A=a['checkpoint_trajectories'][-1][-1]
        B=b['checkpoint_trajectories'][-1][-1]
        if len(b['checkpoint_trajectories'])!=1 or b['time_refinement']!=4*a['time_refinement']:
            raise ArithmeticError('Finer single trial must be twice the previous fine arm')
        if prev is not None and abs(b['starting_scalar_yield']/prev-1)>3e-15:
            raise ArithmeticError('Improved incoming state was not propagated')
        prev=B['scalar_yield']
        if any(x['max_charged']>1 or x['pair_number_error']>3e-13 or abs(x['continuity_residual'])>3e-10
               for x in b['checkpoint_trajectories'][0]):
            raise ArithmeticError('An independent time refinement failed a conservation gate')
        relative=B['scalar_yield']/A['scalar_yield']-1
        if abs(relative)>1e-4:
            raise ArithmeticError(('Time recurrence remains too sensitive',relative))
        rows.append(dict(x=b['x_end'],earlier_fine_yield=A['scalar_yield'],time_refined_yield=B['scalar_yield'],
                         relative_time_change=relative,high_report=b))
    return dict(scope='Complete fine trajectory time comparison on fixed single charge subdivision',stages=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('previous_readout',type=Path)
    p.add_argument('time_refined_segments',nargs=4,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.previous_readout,a.time_refined_segments),indent=2,allow_nan=False))
