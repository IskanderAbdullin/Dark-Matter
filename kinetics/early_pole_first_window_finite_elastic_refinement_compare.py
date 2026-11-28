"""Validate a matched two-projection elastic radial quadrature refinement."""
import argparse
import json
from pathlib import Path
from early_pole_joined_temporal_refinement_reference import conservation


def compare(coarse, fine, refined=None):
    conditions = ('duration', 'steps', 'projection', 'grid', 'aperture_grid',
                  'hard_edge_sets', 'strength', 'same_soft_strength', 'compton_strength',
                  'conversion_strength', 'pair', 'method', 'response_smear',
                  'screening_update', 'elastic_energy_budget')
    out = {}
    for key in ('steps', 'cell-conservative'):
        a, b = coarse[key], fine[key]
        if any(a[k] != b[k] for k in conditions) or a['projection'] != key:
            raise ValueError(f'Mismatched controls for {key}.')
        if any(b[k] != [2, 1, 2] for k in ('unlike_packet_quadrature',
                                            'same_packet_quadrature')):
            raise ValueError(f'Unexpected refined quadrature for {key}.')
        if any(a.get(k, [1, 1, 2]) != [1, 1, 2] for k in
               ('unlike_packet_quadrature', 'same_packet_quadrature')):
            raise ValueError(f'Unexpected baseline quadrature for {key}.')
        if any(a['snapshots'][0][k] != b['snapshots'][0][k] for k in
               ('hard_moments', 'cold_bins', 'instantaneous_source_by_q')):
            raise ValueError(f'Mismatched initial state or cold source for {key}.')
        rates = [p['scalar_gain']/p['duration'] for p in (a, b)]
        out[key] = dict(default=rates[0], radial2=rates[1],
                        change_ppm=1e6*(rates[1]/rates[0]-1),
                        closure=[conservation(p) for p in (a,b)])
    contrast = {c: 100*(1-out['cell-conservative'][c]/out['steps'][c])
                for c in ('default','radial2')}
    if refined is not None:
        for key in ('steps','cell-conservative'):
            a,b=coarse[key],refined[key]
            if any(a[k]!=b[k] for k in conditions) or a['projection']!=key:
                raise ValueError(f'Mismatched radial-four controls for {key}.')
            if any(b[k]!=[4,1,2] for k in ('unlike_packet_quadrature','same_packet_quadrature')):
                raise ValueError(f'Unexpected radial-four graph for {key}.')
            if any(a['snapshots'][0][k]!=b['snapshots'][0][k] for k in
                   ('hard_moments','cold_bins','instantaneous_source_by_q')):
                raise ValueError(f'Mismatched radial-four initial conditions for {key}.')
            rate=b['scalar_gain']/b['duration']
            out[key].update(radial4=rate,
                            change_2_to_4_ppm=1e6*(rate/out[key]['radial2']-1),
                            radial4_closure=conservation(b))
        contrast['radial4']=100*(1-out['cell-conservative']['radial4']/out['steps']['radial4'])
    return dict(cases=out, contrast_percent=contrast,
                scope='Finite first-window sensitivity only; both elastic radial orders changed together.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--default', required=True, nargs=2,
                        help='Default native and continuous raw outputs, in order.')
    parser.add_argument('--radial2', required=True, nargs=2,
                        help='Refined native and continuous raw outputs, in order.')
    parser.add_argument('--radial4',nargs=2,help='Radial-four native and continuous raw outputs.')
    args = parser.parse_args()
    def read(paths):
        return {k:json.loads(Path(p).read_text())
                for k,p in zip(('steps','cell-conservative'), paths, strict=True)}
    print(json.dumps(compare(read(args.default),read(args.radial2),read(args.radial4) if args.radial4 else None),indent=2,allow_nan=False))
