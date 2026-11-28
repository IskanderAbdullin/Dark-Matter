"""Contrast a strict local FRW energy budget at fixed width and a width change at fixed budget."""
import argparse
import json
import math
from pathlib import Path


def summarize(base72, strict72, strict144):
    docs=dict(base72=base72,strict72=strict72,strict144=strict144)
    controls=('projection','grid','duration','expands','final_log_scale','initial_warm_nodes')
    out={}
    for label,case in docs.items():
        active={}
        all_arms=[]
        panels=144 if label=='strict144' else 72
        for expands,shape in enumerate(('fixed','expands')):
            d=case[shape]
            if d['expands']!=bool(expands) or d['panels']!=panels:
                raise ValueError('Wrong profile or FRW background label.')
            if any(d[k]!=base72[shape][k] for k in controls):
                raise ValueError('Mismatched physical controls.')
            if len(d['arms'])!=2:
                raise ValueError('Both warm marker controls are required.')
            for ix,z in enumerate(d['arms']):
                ref=base72[shape]['arms'][ix]
                if any(z[k]!=ref[k] for k in ('initial_count_rates','initial_gap')):
                    raise ValueError('Mismatched initialized warm input.')
                if not math.isclose(z['energy_change']-z['redshift_trace_integral'],
                                    z['redshift_subtracted_defect'],rel_tol=0.,abs_tol=2e-22):
                    raise ValueError('Redshift trace does not reconcile with energy receipt.')
                if z['warm_active']!=bool(ix):
                    raise ValueError('Mismatched warm marker toggle.')
                if not (abs(z['redshift_subtracted_defect'])<3e-12
                        and abs(z['entropy_defect'])<3e-11 and abs(z['hard_count_defect'])<3e-12):
                    raise ValueError('Terminal conservation gate failed.')
                if label.startswith('strict') and (z['elastic_energy_budget']!=5e-15 or
                        z['peak_abs_elastic_energy_defect']>=z['elastic_energy_budget']):
                    raise ValueError('Local budget gate failed.')
            active[shape]=d['arms'][1]
            all_arms.extend(d['arms'])
            if any(not math.isfinite(active[shape][k]) or active[shape][k]<=0.
                   for k in ('cold_gain','warm_gain')):
                raise ValueError('Invalid active-marker growth for paired difference.')
        same=('warm_gain','cold_gain')
        redshift={k:1e6*(active['expands'][k]/active['fixed'][k]-1) for k in same}
        out[label]=dict(redshift_ppm=redshift,active=active,
                        max_abs_trace_subtracted_energy_defect_4_arms=max(
                            abs(x['redshift_subtracted_defect']) for x in all_arms),
                        max_abs_entropy_defect_4_arms=max(
                            abs(x['entropy_defect']) for x in all_arms))
    def difference(a,b):
        return {k:{field:1e6*(out[b]['active'][k][field]/out[a]['active'][k][field]-1)
                   for field in ('cold_gain','warm_gain')} for k in ('fixed','expands')}
    return dict(cases=out,change_budget_at_72_ppm=difference('base72','strict72'),
                change_width_at_fixed_budget_ppm=difference('strict72','strict144'),
                scope='Finite profile differences, not a two-axis interaction or convergence bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for profile in ('base72','strict72','strict144'):
        p.add_argument('--'+profile,nargs=2,required=True,help='Fixed, then expanding output.')
    a=p.parse_args()
    def read(files):
        return {k:json.loads(Path(f).read_text()) for k,f in zip(('fixed','expands'),files,strict=True)}
    print(json.dumps(summarize(read(a.base72),read(a.strict72),read(a.strict144)),indent=2,allow_nan=False))
