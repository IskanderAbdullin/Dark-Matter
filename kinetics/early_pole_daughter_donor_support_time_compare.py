"""Compare two admissible equal-bare-moment retunes on one fixed held chart."""
import argparse
import json
from pathlib import Path

import numpy as np


KEYS=('native','local-daughter','retuned-local-daughter','retuned-local-daughter-no-neighbors')


def compare(data,static):
    native=data[0]
    result={}
    h={}
    for v,k in zip(data,KEYS,strict=True):
        if v.get('prepared_charged_swap','native')!=k:
            raise ValueError('Unexpected preparation label.')
        for f in ('duration','projection','rtol','base_atol','history_exposure','report_panels'):
            if v[f]!=native[f]:
                raise ValueError('Numerical controls differ.')
        if v.get('elastic_channels',True) is not True or v.get('selected_elastic_channel','both')!='both':
            raise ValueError('Require both held elastic channels.')
        if not (abs(v['energy_defect'])<3e-12 and abs(v['entropy_defect'])<3e-11
                and abs(v['hard_number_defect'])<3e-12):
            raise ValueError('Unexpected terminal defect.')
        c=np.asarray(v['cold_count_history'])
        if c.shape!=(native['report_panels']+1,) or c[0]!=native['cold_count_history'][0]:
            raise ValueError('Different recorded grid or cold seed.')
        h[k]=c-c[0]
    for v,label in zip(data[2:],('all_outside','no_neighbors'),strict=True):
        setup=v['complement_retune']
        actual=static['cases'][label]
        for k in ('fixed_cells','complement_cells','number_and_bare_energy_defect','affine_logit_multipliers'):
            if setup[k]!=actual[k]:
                raise ValueError('The static and evolved inputs differ.')
        if max(map(abs,setup['number_and_bare_energy_defect']))>=5e-14:
            raise ValueError('Bare daughter moment gate failed.')
        result[label]=setup
    if result['all_outside']['fixed_cells']!=result['no_neighbors']['fixed_cells']:
        raise ValueError('Compensation changed the aperture-local intervention.')
    outside=set(result['all_outside']['complement_cells'])
    excluded=outside-set(result['no_neighbors']['complement_cells'])
    if len(excluded)!=2:
        raise ValueError('Expected removal of both immediate outside cells.')
    t=native['duration']
    a,b=(h[k] for k in KEYS[2:])
    return dict(mean_rate={k:float(v[-1]/t) for k,v in h.items()},
                local_mean_effect={k:float((v[-1]-h['native'][-1])/t) for k,v in h.items()},
                terminal_mean_change_no_neighbors_vs_all=float((b[-1]-a[-1])/t),
                cumulative_change_no_neighbors_vs_all=(b-a).tolist(),
                exposure=native['history_exposure'],retune_settings=result,
                excluded_donor_indices=sorted(excluded),
                scope='Same two bare daughter moments and local data; different permitted outside shapes only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','local','all_outside','no_neighbors','static'):
        p.add_argument(k)
    a=p.parse_args()
    values=list(vars(a).values())
    print(json.dumps(compare([json.loads(Path(f).read_text()) for f in values[:4]],
                             json.loads(Path(values[4]).read_text())),
                     indent=2,allow_nan=False))
