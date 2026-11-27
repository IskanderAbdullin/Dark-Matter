"""Compact diagnostics across paired finite-volume solver controls."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np


def go(root,archive=None):
    roots=['small','fine','small-direct','fine-direct','vacuum-direct','only-pauli','finite-daughter']
    rows={}
    for s in roots:
        path=root/(s+'.json')
        raw=path.read_bytes()
        q=json.loads(raw)
        f=q['checkpoints'][-1]
        assert q['empty_daughter_row_sum_control']<1e-10
        assert f['max_fermion_occupation']<1.
        assert abs(sum(q['final_scalar_n_per_visible_entropy'])/f['Y']-1)<3e-14
        rows[s]=dict(final=f,worst_number_defect=q['largest_sampled_pair_number_defect'],
                     scalar_bins=q['scalar_bins'],fermion_bins=q['fermion_bins'],solver_calls=q['solver_calls'])
        if archive is not None:
            name=f'early-pole-paired-hartree-{s}-2026-11-26.json.gz'
            (archive/name).write_bytes(gzip.compress(raw,mtime=0))
    def y(s):
        return rows[s]['final']['Y']
    c=abs(y('small-direct')/y('fine-direct')-1.)
    inverse=y('small-direct')/y('only-pauli')
    assert c<.018
    assert max(abs(y(s)/y(s+'-direct')-1.) for s in ('small','fine'))<4e-6
    assert max(rows[s]['worst_number_defect'] for s in roots if s not in ('small','fine'))<2e-7
    return dict(scope="Common-event streaming-only held-parent comparison. No reaction including "
                      "energy-resolved parent evolution, actual reheating, additional scattering or field dynamics.",
                paired_vs_empty_daughter_previous_fine_yield=y('fine-direct')/.0008022372208537028,
                inverse_enabled_over_Pauli_only_on_fixed_coarse_grids=inverse,
                combined_momentum_grid_relative_change=c,runs=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run_directory',type=Path)
    p.add_argument('--archive',type=Path)
    a=p.parse_args()
    print(json.dumps(go(a.run_directory,a.archive),indent=2,allow_nan=False))
