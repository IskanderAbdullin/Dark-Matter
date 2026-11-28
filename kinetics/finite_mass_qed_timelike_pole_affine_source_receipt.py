"""Validate matched fixed-angle transverse and Coulomb coefficient meshes."""
import argparse
import json
from pathlib import Path

import numpy as np


def run(root):
    outcome=[]
    expected={('T',.0015),('T',.00256),('L',.0012),('L',.00256)}
    for g in ('fine','third'):
        rec=[]
        for mesh in (0,1):
            p=root/f'timelike-pole-{g}-mesh{mesh}.jsonl'
            with p.open() as stream:
                rows=[json.loads(s) for s in stream]
            normalized={(r['polarization'],r['k']):r for r in rows}
            if len(rows)!=4 or set(normalized)!=expected or any(r['mesh']!=([3,6] if mesh==0 else [6,12])
                      or r['grid']!=g or r['c']!=0. for r in rows):
                raise ArithmeticError('Mismatched archived pole coefficients')
            rec.append(normalized)
        for key in sorted(expected):
            a,b=[r[key] for r in rec]
            v=np.array(a['raw_kernel'])
            w=np.array(b['raw_kernel'])
            z=float(max(abs(v/w-1)))
            if any(w<=0.) or z>8e-7 or min(a['pole_population_neutral'],b['pole_population_neutral'])<=0.:
                raise ArithmeticError(('Fixed aperture mesh',g,key,z))
            b['joint_phase_and_azimuth_max_mesh_relative_defect']=z
            b['neutral_population_max_mesh_relative_defect']=abs(a['pole_population_neutral']/
                                                                b['pole_population_neutral']-1)
            outcome.append(b)
    return dict(scope='Fixed-k,c conditional pole absorption/inverse kernel only, bare vertices with retarded frozen residues.',
                data=outcome)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
