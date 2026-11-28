"""Exercise accepted-step restoration in the published frozen RHS (in-process audit)."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

import numpy as np

import early_pole_joined_simultaneous_chart_reference as base
from accepted_dop853_checkpoint_reference import evolve,dense_matrix

seed='results/early-pole-first-window-hard1-cold32-paired-initializer-2026-11-27.json.gz'
saved=base.solve_ivp

def replay(f,tspan,y,method,rtol,atol,max_step,t_eval):
    if method!='DOP853':raise ValueError(method)
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'adaptive.json'
        kwargs=(tspan,y,t_eval,max_step,rtol,atol,'frozen-chart-H1',path)
        if evolve(f,*kwargs,stop=1) is not None:
            raise RuntimeError('Fixture did not suspend at accepted step 1')
        b=evolve(f,*kwargs,start=path)
    return SimpleNamespace(y=dense_matrix(b),success=True, message='validated in-process restore')

args=(seed,'steps',1.25e-10,2,3e-9,2e-14)
x=base.run(*args,record_mode_history=True)
try:
    base.solve_ivp=replay
    y=base.run(*args,record_mode_history=True)
finally:
    base.solve_ivp=saved
if x!=y:
    k=[z for z in x if x[z]!=y.get(z)]
    raise ArithmeticError(('Adaptive frozen replay not bitwise',k))
print(json.dumps(dict(scope='H1 frozen both-elastic small-step in-process RHS checkpoint only',same_raw_result=True,
                      original_calls=x['calls'],duration=x['duration'],kept_panels=x['report_panels'],
                      explicit_energy=x['energy_defect'],explicit_entropy=x['entropy_defect']),indent=2))
