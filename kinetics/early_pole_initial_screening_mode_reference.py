"""Compare held, static-only and fully refreshed screening on the same early time grid."""
import argparse
import json

import numpy as np

from early_pole_joined_temporal_refinement_reference import conservation, select


def run(runs):
    ref=runs['frozen']
    guards=('grid','projection','pair','method','strength','compton_strength','same_soft_strength',
            'conversion_strength','duration','steps','hard_edge_sets','aperture_grid','response_smear',
            'elastic_energy_budget')
    table={}
    comparison={}
    for k,x in runs.items():
        if k!=x['screening_update'] or any(x[x0]!=ref[x0] for x0 in guards):
            raise ValueError('The screening mode must be the only changed numerical control.')
        for j in ('hard_moments','cold_bins','instantaneous_source_by_q'):
            np.testing.assert_array_equal(ref['snapshots'][0][j],x['snapshots'][0][j])
        z=np.asarray(x['state_history'])
        q=np.asarray(ref['state_history'])
        if z.size!=q.size or z.size!=x['steps']+1 or not np.all(np.isfinite(z)) or z[0]!=q[0]:
            raise ValueError('Use complete same-grid trajectories.')
        if z[-1]-z[0]!=x['scalar_gain']:
            raise ValueError('Wrong state history gain.')
        row=x['response_history']
        if k=='frozen':
            assert not row
        else:
            assert len(row)==x['steps']
        table[k]=dict(mean_rate=x['scalar_gain']/x['duration'],
                      gain_relative_to_held=x['scalar_gain']/ref['scalar_gain']-1.,
                      supremum_absolute_count_difference_from_held=float(np.max(abs(z-q))),
                      terminal_source=float(sum(x['snapshots'][-1]['instantaneous_source_by_q'])),
                      maximum_absolute_debye_squared_fractional_change=max(
                          (abs(r['fractional_mD2_change']) for r in row),default=0.),
                      maximum_absolute_soft_packet_ratio_shift=max(
                          (abs(v-1.) for r in row for v in r['soft_packet_ratio_quantiles']),default=0.),
                      conservation=conservation(x))
        comparison[k]=(z-q).tolist()
    return dict(projection=ref['projection'],duration=ref['duration'],steps=ref['steps'],
                compared=table,absolute_count_differences_from_held=comparison,
                full_minus_debye_gain_fraction=(runs['full']['scalar_gain']-runs['debye']['scalar_gain'])/ref['scalar_gain'],
                scope='Matching finite early seeded histories. Same-flavor polarization held on all arms.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    for a in ('frozen','debye','full'):
        p.add_argument('--'+a,required=True)
    args=vars(p.parse_args())
    pp=args.pop('projection')
    print(json.dumps(run({k:select(v,pp) for k,v in args.items()}),indent=2,allow_nan=False))
