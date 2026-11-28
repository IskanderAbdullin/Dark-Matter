"""Compare the local prepared chart on original eight and nested thirty-two cold bins."""
import argparse
import json
from pathlib import Path

import numpy as np


def validate(pair, n):
    a,b=pair
    if a.get('prepared_charged_swap','native')!='native' or b.get('prepared_charged_swap')!='local-daughter':
        raise ValueError('Expected native and exact aperture-local prepared profiles.')
    controls=('duration','rtol','base_atol','projection','report_panels','history_exposure')
    if any(a[k]!=b[k] for k in controls):
        raise ValueError('Unexpected chart controls.')
    for x in (a,b):
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError('Original closure gate failed.')
        if (x.get('aperture_subdivision',1)!=1 or x.get('elastic_channels',True) is not True
                or x.get('selected_elastic_channel','both')!='both'):
            raise ValueError('Only the original hard mesh and full held elastic terms are comparable.')
    E=np.geomspace(1e-6,1e-5,n+1)
    q=np.cbrt((E[1:]**3+E[:-1]**3)/2)
    rows=[]
    for x in (a,b):
        c=np.asarray(x['cold_mode_count_history'])
        if c.shape!=(a['report_panels']+1,n):
            raise ValueError('Wrong cold-bin width.')
        np.testing.assert_allclose(x['cold_momentum_nodes_over_T'],q,rtol=2e-15,atol=0.)
        np.testing.assert_allclose(c.sum(axis=1),x['cold_count_history'],rtol=0.,atol=1e-16)
        rows.append(c)
    np.testing.assert_array_equal(rows[0][0],rows[1][0])
    return rows[0][0],rows[1]-rows[0],np.asarray(rows),E,controls


def compare(original, refined):
    r0,d0,c0,e0,controls=validate(original,8)
    r1,d1,c1,e1,_=validate(refined,32)
    if any(original[0][k]!=refined[0][k] for k in controls):
        raise ValueError('Mismatch beyond cold grid refinement.')
    np.testing.assert_allclose(e1[::4],e0,rtol=2e-15,atol=0.)
    np.testing.assert_allclose(r1.reshape(8,4).sum(axis=1),r0,rtol=0.,atol=4e-17)
    dgroup=d1.reshape(len(d1),8,4).sum(axis=2)
    exp=original[0]['duration']
    a0=(c0[:,-1]-c0[:,0])/exp
    a1=((c1[:,-1]-c1[:,0])/exp).reshape(2,8,4).sum(axis=2)
    return dict(absolute_means_native_then_local=[a0.sum(axis=1).tolist(),a1.sum(axis=1).tolist()],
                grouped_absolute_mean_change_native_then_local=(a1-a0).tolist(),
                max_initial_group_seed_difference=float(np.max(abs(r1.reshape(8,4).sum(axis=1)-r0))),
                terminal_original_effect=(d0[-1]/exp).tolist(),
                terminal_refined_effect=(d1[-1]/exp).tolist(),
                terminal_refined_aggregated_effect=(dgroup[-1]/exp).tolist(),
                terminal_aggregate_gain_change_by_mode=((dgroup[-1]-d0[-1])/exp).tolist(),
                total_terminal_mean_effect=[float(d[-1].sum()/exp) for d in (d0,d1)],
                original_final_negative_indices=np.flatnonzero(d0[-1]<0).tolist(),
                refined_final_negative_indices=np.flatnonzero(d1[-1]<0).tolist(),
                exposure=original[0]['history_exposure'],
                original_cumulative_effect=d0.tolist(),refined_aggregated_cumulative_effect=dgroup.tolist(),
                scope='Nested cold-bin refinement for a fixed held hard graph, not a continuum bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('original','refined'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native, then aperture-local histories.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in paths]
                                 for paths in vars(a).values())),indent=2,allow_nan=False))
