"""Refine only the held elastic stage within an unchanged outer split schedule."""
import argparse
import json
from math import fsum
from unittest.mock import patch

import early_pole_cold_screened_split_reference as split
from early_pole_prepared_hard_swap_trajectory_reference import run as prepared


def run(bundle, case, macro_steps, elastic_substeps):
    if macro_steps<1 or elastic_substeps<1:
        raise ValueError('Positive outer and elastic step counts are required.')
    original=split.step

    def subcycle(model, y, dt, method='backward', energy_budget=None):
        if energy_budget!=5e-15 or method!='backward':
            raise ValueError('Use the matched local elastic energy gate.')
        ledgers=[]
        for _ in range(elastic_substeps):
            y,r=original(model,y,dt/elastic_substeps,method,energy_budget=energy_budget)
            if abs(r['rho_drift'])>=energy_budget:
                raise ValueError('Recorded small elastic step exceeded the local budget.')
            ledgers.append(r)
        return y,dict(h=dt,elastic_substeps=ledgers,
                      rho_drift=fsum(x['rho_drift'] for x in ledgers),
                      entropy_gain=fsum(x['entropy_gain'] for x in ledgers))

    with patch.object(split,'step',subcycle):
        z=prepared(bundle,case,macro_steps,hold_screening=True)
    z['elastic_substep_energy_budget']=z.pop('elastic_energy_budget')
    z['elastic_substeps_per_macro_step']=elastic_substeps
    return z


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle')
    p.add_argument('--case',choices=('native','local-daughter'),required=True)
    p.add_argument('--macro-steps',default=16,type=int)
    p.add_argument('--elastic-substeps',default=4,type=int)
    q=p.parse_args()
    print(json.dumps(run(**vars(q)),indent=2,allow_nan=False))
