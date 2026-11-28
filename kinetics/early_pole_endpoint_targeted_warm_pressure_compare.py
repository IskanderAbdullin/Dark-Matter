"""Compare the evolved targeted warm chart with independently evolved uniform charts."""
import argparse
import json
from math import fsum
from pathlib import Path
import numpy as np
from early_pole_endpoint_gauss_warm_comparison import gauss,chart_for
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_pressure_cells import partition


def run(base,target,uniforms):
    x=base['x_end']
    d,y=unpack(base,False,x)
    m=d.state(x,y)[0]
    if chart_for(target)[0]!='targeted-2xG4' or [chart_for(q)[0] for q in uniforms] not in ([4],[4,5]):
        raise ArithmeticError('Expected targeted and independently evolved matched G4 and optional G5')
    tv,env,C=gauss(target,d,m,x)
    selected=np.isin(d.ds.ed[:-1],target['targeted_rule']['target_q_lower_edges'])
    if np.count_nonzero(selected)!=2:
        raise ArithmeticError('Missing targeted parent cells')
    out=[]
    for card in uniforms:
        uv,renv,S=gauss(card,d,m,x)
        r=partition(S,C)/uv[1]
        ref=uv[1]
        if abs(fsum(r[:,0])-(tv[1]/ref-1))>4e-15:
            raise ArithmeticError('Pressure difference did not close on common parents')
        p,q=env['pressure_total_min_max']
        l,h=renv['pressure_total_min_max']
        out.append(dict(uniform_order=card['warm_gauss_order'],target_over_uniform_minus_one=(tv/uv-1).tolist(),
                        pressure_ratio_interval=[p/h-1,q/l-1],fixed_reference_uniform_pressure=ref,
                        all_pressure_differences=[fsum(k) for k in r.T],
                        targeted_parents_only_pressure_differences=[fsum(k) for k in r[selected].T],
                        untouched_parents_pressure_differences=[fsum(k) for k in r[~selected].T]))
    return dict(fixed_x=x,fixed_mass=m,targeted_moments=tv.tolist(),comparisons=out,
                moment_columns=['number','pressure','warm_pressure'],
                pressure_difference_columns=['total','changed_cell_count','changed_pressure_per_particle'],
                scope='Separate conditional finite evolutions at a common measurement mass; no claim targeted improves physical accuracy.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('base',type=Path)
    p.add_argument('target',type=Path)
    p.add_argument('uniforms',type=Path,nargs='+')
    a=p.parse_args()
    print(json.dumps(run(read(a.base),read(a.target),[read(x) for x in a.uniforms]),indent=2,allow_nan=False))
