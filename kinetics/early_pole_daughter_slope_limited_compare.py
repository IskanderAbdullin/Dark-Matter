"""Audit the sampled slope-cap range and resolved signs without extrapolating."""
import argparse
import json
from pathlib import Path


def compare(runs,thresholds,dual_checks):
    detail=[]
    maxdefect=0.
    specs=((4,64,True,False),(4,128,True,False),(4,64,False,False),
           (2,256,True,False),(4,128,True,True),(1,512,True,False))
    if len(runs)!=len(specs):
        raise ValueError('Required six resolution and constraint controls.')
    for i,(mode,n,gauge,exterior) in enumerate(specs):
        if len(runs[i])!=mode:
            raise ValueError('Wrong number of selected hard/cold fixtures.')
        for v in runs[i]:
            if (v['fixed_gauge_moment']!=gauge or v['subdivision']!=n
                    or ('fixed_exterior_occupations' in v)!=exterior):
                raise ValueError('Mismatched slope or constraint profile.')
            slopes=[]
            failures=[]
            positive=[]
            outer_negative=[]
            for row in v['cases']:
                L=row['slope_cap']
                if L<=0. or (slopes and L<=slopes[-1]):
                    raise ValueError('Require increasing positive trial slopes.')
                slopes.append(L)
                pair=row['paired_witness_range']
                if pair is not None:
                    for arm in ('native','local'):
                        for tag in ('minimum','maximum'):
                            p=row[arm][tag]
                            if p['verified_max_slope']>L:
                                raise ValueError('A trial occupation exceeds the claimed slope.')
                            maxdefect=max(maxdefect,p['normalized_moment_defect'])
                    if pair[1]>0.:
                        positive.append(L)
                elif any(row[k]['status']=='numerical_witness_repair_failed' for k in ('native','local')):
                    failures.append(L)
                o=row.get('numerical_interpolation_outer_range')
                if o is not None and o[1]<0.:
                    outer_negative.append(L)
            detail.append(dict(hard_depth=v['hard_depth'],cold_modes=v['cold_modes'],subdivision=n,
                               fixed_gauge_moment=gauge,fixed_exterior=exterior,first_tested_positive_slope=min(positive) if positive else None,
                               largest_tested_negative_outer_slope=max(outer_negative) if outer_negative else None,
                               explicit_numerical_repair_failures=failures))
    if len(thresholds)!=3 or any(len(z)!=4 for z in thresholds):
        raise ValueError('Need both real slope objectives for every constraint set.')
    for cases,gauge,exterior in zip(thresholds,(True,True,False),(False,True,False),strict=True):
        for a in cases:
            if a['fixed_gauge_moment']!=gauge or a['match_exterior']!=exterior:
                raise ValueError('The reported constraints differ from the matched test.')
            g=a['minimum_feasibility']
            q=a['minimum_positive']
            if (q['verified_paired_rate']<=0. or
                    g['numerical_finite_lower_slope']>q['verified_witness_slope']+1e-8 or
                    max(a['analytic_average_based_feasibility_lower_slope'])>g['verified_witness_slope']+1e-8):
                raise ValueError('A finite threshold does not satisfy independent gates.')
    report_outer={}
    for batch in dual_checks:
        for v in batch:
            key=f"{v['hard_depth']}x{v['cold_modes']}|{v['subdivision']}"
            selected=[r for r in v['cases'] if r['numerical_interpolation_outer_range'] is not None]
            report_outer[key]={str(r['slope_cap']):r['numerical_interpolation_outer_range'] for r in selected}
    return dict(scope='Slope in occupation units per unit of p/T. No inference outside the frozen model.',
                checked=detail,largest_normalized_moment_defect=maxdefect,
                weak_dual_outer_checks=report_outer,thresholds=thresholds,raw=runs,raw_dual_checks=dual_checks)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs=6)
    p.add_argument('--thresholds',nargs=3,required=True)
    p.add_argument('--dual-checks',nargs='+',required=True)
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in getattr(a,k)] for k in ('inputs','thresholds','dual_checks'))),
                     indent=2,allow_nan=False))
