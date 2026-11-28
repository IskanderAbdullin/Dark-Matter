"""Compare matched held-polarization warm FRW runs under two numerical profiles."""
import argparse
import json
from pathlib import Path


def compare(old, strict):
    cases={}
    for key, expands in (('fixed',False),('expands',True)):
        a,b=old[key],strict[key]
        common=('projection','grid','duration','expands','final_log_scale','initial_warm_nodes')
        if any(a[k]!=b[k] for k in common) or b['expands']!=expands:
            raise ValueError(f'Mismatched FRW controls: {key}.')
        if b['panels']!=2*a['panels'] or b.get('collision_clock','radiation')!='radiation':
            raise ValueError(f'Mismatched FRW refinement or clock: {key}.')
        for idx,(x,y) in enumerate(zip(a['arms'],b['arms'],strict=True)):
            if x['warm_active']!=bool(idx) or y['warm_active']!=bool(idx):
                raise ValueError(f'Mismatched FRW marker toggle: {key}.')
            if abs(y['redshift_subtracted_defect'])>=3e-12 or abs(y['entropy_defect'])>=3e-11:
                raise ValueError(f'Terminal closure rejected: {key}, arm {idx}.')
            if not 0<=y['peak_abs_elastic_energy_defect']<y['elastic_energy_budget']==5e-15:
                raise ValueError(f'Elastic step closure rejected: {key}, arm {idx}.')
        fields=('cold_gain','warm_gain','hard_count_defect','redshift_trace_integral',
                'redshift_subtracted_defect','entropy_defect','elastic_splits',
                'peak_abs_elastic_energy_defect','sum_absolute_elastic_energy_defects')
        cases[key]={profile: [{k:z[k] for k in fields if k in z} for z in obj['arms']]
                    for profile,obj in (('base',a),('strict',b))}
    effects={}
    for p,docs in (('base',old),('strict',strict)):
        a,b=(docs[k]['arms'][1] for k in ('fixed','expands'))
        effects[p]={k:1e6*(b[k]/a[k]-1) for k in ('cold_gain','warm_gain')}
    return dict(cases=cases, expanding_minus_fixed_ppm=effects,
                profile_panels=[old['fixed']['panels'],strict['fixed']['panels']],
                scope='Simultaneous macro step and local energy budget change under a prescribed radiation clock.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',nargs=2,required=True,help='Fixed and expanding runs in order.')
    p.add_argument('--strict',nargs=2,required=True,help='Fixed and expanding runs in order.')
    a=p.parse_args()
    def read(paths):
        return {k:json.loads(Path(f).read_text()) for k,f in zip(('fixed','expands'),paths,strict=True)}
    print(json.dumps(compare(read(a.base),read(a.strict)),indent=2,allow_nan=False))
