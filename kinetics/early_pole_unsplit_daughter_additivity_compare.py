"""Compare a complete held-graph daughter shape factorial and temporal nonadditivity."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(raw_inputs, histories):
    labels=('native','local-daughter','outside-daughter','daughter')
    if len(raw_inputs)!=4 or len(histories)!=4:
        raise ValueError('Require four inputs and trajectories.')
    ys=[np.asarray(a['snapshots'][0]['hard_moments']) for a in raw_inputs]
    a,b,c,d=ys
    m1=np.any(b!=a,axis=1)
    m2=np.any(c!=a,axis=1)
    if np.any(m1&m2) or int(m1.sum())!=2 or int(m2.sum())!=18:
        raise ValueError('Does not have the specified disjoint hard input supports.')
    np.testing.assert_array_equal(np.where(m1[:,None],b,c),d)
    for i,x in enumerate(raw_inputs):
        if i and x.get('prepared_charged_swap')!=labels[i]:
            raise ValueError('Wrong preparation input.')
        for k in ('projection','grid','aperture_grid','hard_edge_sets'):
            if x[k]!=raw_inputs[0][k]:
                raise ValueError(f'Mismatched setup: {k}.')
        np.testing.assert_array_equal(x['snapshots'][0]['cold_bins'],
                                      raw_inputs[0]['snapshots'][0]['cold_bins'])
    history=[]
    for i,x in enumerate(histories):
        if i and x.get('prepared_charged_swap')!=labels[i]:
            raise ValueError('Wrong trajectory label.')
        for k in ('projection','duration','rtol','base_atol','report_panels','history_exposure'):
            if x[k]!=histories[0][k]:
                raise ValueError(f'Mismatched solver or time grid: {k}.')
        if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                and abs(x['hard_number_defect'])<3e-12):
            raise ValueError('Unaccepted unsplit closure.')
        h=np.asarray(x['cold_count_history'])
        if h.shape!=(x['report_panels']+1,) or h[0]!=histories[0]['cold_count_history'][0]:
            raise ValueError('Unexpected count history shape or initial count.')
        history.append(h-h[0])
    a,b,c,d=history
    cross=(d-a)-(b-a)-(c-a)
    duration=histories[0]['duration']
    return dict(mean_rate={k:float(x[-1]/duration) for k,x in zip(labels,history,strict=True)},
                whole_effect=float((d[-1]-a[-1])/duration),
                local_effect=float((b[-1]-a[-1])/duration),
                outside_effect=float((c[-1]-a[-1])/duration),
                interaction_rate=float(cross[-1]/duration),interaction_count=cross.tolist(),
                exposure=histories[0]['history_exposure'],
                scope='Finite held-graph prepared-input interaction on one hard mesh only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw-inputs',nargs=4,required=True)
    p.add_argument('--histories',nargs=4,required=True)
    a=p.parse_args()
    def read(paths):
        return [json.loads(Path(f).read_text()) for f in paths]
    print(json.dumps(compare(read(a.raw_inputs),read(a.histories)),indent=2,allow_nan=False))
