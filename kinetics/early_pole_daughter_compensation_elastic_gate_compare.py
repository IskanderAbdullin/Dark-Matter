"""Measure the prepared complement response with and without the elastic vector field."""
import argparse
import json
from pathlib import Path

import numpy as np


LABELS=('native','local-daughter','local-complement-only','retuned-local-daughter')
CONTROLS=('projection','duration','report_panels','rtol','base_atol','history_exposure')


def compare(joined, nonelastic):
    base=joined[0]
    receipt={}
    for label,group in (('joined',joined),('nonelastic',nonelastic)):
        paths={}
        for key,x,ref in zip(LABELS,group,joined,strict=True):
            if x.get('prepared_charged_swap','native')!=key:
                raise ValueError('Wrong prepared input order.')
            if any(x[k]!=base[k] for k in CONTROLS):
                raise ValueError('Mismatched numerical controls.')
            if x.get('elastic_channels',True)!=(label=='joined'):
                raise ValueError('Mismatched elastic switch.')
            if x.get('complement_retune')!=ref.get('complement_retune'):
                raise ValueError('Complement preparation changed with channel switch.')
            if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                    and abs(x['hard_number_defect'])<3e-12):
                raise ValueError('Original terminal gate failed.')
            z=np.asarray(x['cold_count_history'])
            if z.shape!=(base['report_panels']+1,) or z[0]!=base['cold_count_history'][0]:
                raise ValueError('Mismatched cold seed or exposure width.')
            paths[key]=z-z[0]
        dur=base['duration']
        n,a,b,c=(paths[k] for k in LABELS)
        receipt[label]=dict(final_rate={k:float(v[-1]/dur) for k,v in paths.items()},
                        local_effect=float((a[-1]-n[-1])/dur),
                        complement_effect=float((b[-1]-n[-1])/dur),
                        combined_effect=float((c[-1]-n[-1])/dur),
                        finite_interaction=float(((c-n)-(a-n)-(b-n))[-1]/dur),
                        compensation_on_local_arm_effect=float((c[-1]-a[-1])/dur),
                        cumulative_complement_effect=(b-n).tolist(),
                        cumulative_added_compensation_effect=(c-a).tolist())
    a,b=(receipt[k] for k in ('joined','nonelastic'))
    return dict(cases=receipt,exposure=base['history_exposure'],
                complement_only_without_over_with=b['complement_effect']/a['complement_effect'],
                added_compensation_without_over_with=(
                    b['compensation_on_local_arm_effect']/a['compensation_on_local_arm_effect']),
                scope='Finite conditional switch at held dielectric tables; not an additive physical attribution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('joined','nonelastic'):
        p.add_argument('--'+k,required=True,nargs=4,help='Native, local, compensation only, both.')
    a=p.parse_args()
    def read(paths):
        return [json.loads(Path(x).read_text()) for x in paths]
    print(json.dumps(compare(read(a.joined),read(a.nonelastic)),indent=2,allow_nan=False))
