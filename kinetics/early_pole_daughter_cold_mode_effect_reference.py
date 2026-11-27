"""Resolve the prepared daughter count effect mode by mode at saved split snapshots."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation


def project(native, swapped):
    controls=('projection','grid','duration','steps','aperture_grid','hard_edge_sets',
              'strength','same_soft_strength','compton_strength','conversion_strength',
              'screening_update','method','pair','response_smear','elastic_energy_budget')
    if any(native[k]!=swapped[k] for k in controls):
        raise ValueError('Mismatched split controls.')
    if swapped['prepared_charged_swap'] not in ('local-daughter','outside-daughter'):
        raise ValueError('Require a local or outside-only daughter swap.')
    a,b=native['snapshots'],swapped['snapshots']
    exp=np.array([x['exposure'] for x in a])
    np.testing.assert_array_equal(exp,[x['exposure'] for x in b])
    y,z=(np.array([x['cold_bins'] for x in data]) for data in (a,b))
    np.testing.assert_array_equal(y[0],z[0])
    gain=(z[1:]-y[1:])/exp[1:,None]
    total=gain.sum(axis=1)
    absolute=np.abs(gain).sum(axis=1)
    cancel=np.divide(absolute-np.abs(total),absolute,
                     out=np.zeros_like(total),where=absolute>0.)
    initial=(np.asarray(b[0]['instantaneous_source_by_q'])
             -np.asarray(a[0]['instantaneous_source_by_q']))
    ref=swapped['scalar_gain']/swapped['duration']-native['scalar_gain']/native['duration']
    np.testing.assert_allclose(total[-1],ref,atol=8e-11,rtol=0.)
    return dict(steps=native['steps'],screening_update=native['screening_update'],
                preparation=swapped['prepared_charged_swap'],exposure=exp[1:].tolist(),
                mean_mode_effects=gain.tolist(),mean_mode_sum=total.tolist(),
                mean_mode_l1=absolute.tolist(),
                cancellation_fraction=cancel.tolist(),
                initial_mode_shifts=initial.tolist(),conservation=[conservation(native),conservation(swapped)],
                scope='Saved cold-bin count differences only; no derivative or mode convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native')
    p.add_argument('swap')
    a=p.parse_args()
    print(json.dumps(project(*(json.loads(Path(p).read_text()) for p in vars(a).values())),
                     indent=2,allow_nan=False))
