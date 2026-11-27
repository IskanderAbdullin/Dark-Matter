"""Compare the actual half-aperture trajectory with an endpoint interpolation."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def compare(depth1,depth4,static):
    c={}
    th=depth1[0]['duration']
    controls=('projection','duration','report_panels','rtol','base_atol','history_exposure')
    x=depth1[0]
    for depth,raw in ((1,depth1),(4,depth4)):
        seed=x['cold_count_history'][0]
        paths=[]
        for v,k in zip(raw,('native','local-daughter-fractional','local-daughter'),strict=True):
            if (any(v[a]!=x[a] for a in controls) or v.get('aperture_subdivision',1)!=depth
                    or v.get('prepared_charged_swap','native')!=k or not v.get('elastic_channels',True)
                    or v.get('selected_elastic_channel','both')!='both'):
                raise ValueError('Mismatched geometry or kinetics beyond replacement fraction.')
            if not (abs(v['energy_defect'])<3e-12 and abs(v['entropy_defect'])<3e-11 and abs(v['hard_number_defect'])<3e-12):
                raise ValueError('An earlier solver gate was not met.')
            n=np.asarray(v['cold_count_history'])
            if n.shape!=(x['report_panels']+1,) or n[0]!=seed or not np.all(np.isfinite(n)&(n>0.)):
                raise ValueError('Unexpected count history.')
            paths.append(n)
        old,middle,young=paths
        t=np.asarray(x['history_exposure'])
        ref=next(o for o in static if o['depth']==depth and o['cold_modes']==8)
        half=raw[1]
        h=half['local_archive_support']
        if h['replacement_fraction']!=.5 or h['indices']!=ref['indices']:
            raise ValueError('Wrong fractional prepared state.')
        np.testing.assert_allclose(h['momentum_interval'],[1.4138888888888888,1.4666666666666666],rtol=0.,atol=2e-14)
        if not np.isclose(h['initialized_cold_rate'],sum(ref['local_rate_by_mode'][4]),rtol=0.,atol=1e-16):
            raise ValueError('The fractional initializer does not match the static check.')
        v=np.asarray(half['cold_mode_count_history'])
        if v.shape!=(len(t),8) or not np.all(np.isfinite(v)&(v>0.)):
            raise ValueError('Missing accepted half-case cold-mode history.')
        if np.max(abs(v.sum(axis=1)-middle)/np.spacing(middle))>4.:
            raise ValueError('The recorded half-case modes do not sum to the scalar number.')
        tip=middle-old
        end=young-old
        y=tip-end/2
        fn=lambda a:[t[i:i+2].tolist() for i in np.flatnonzero((a[:-1]<0.)&(a[1:]>=0.))]
        c[str(depth)]=dict(half_mean=float(tip[-1]/th),full_mean=float(end[-1]/th),
                           half_relative_to_endpoint_linear_mean=float(y[-1]/(end[-1]/2)),
                           half_mismatched_from_endpoint_linear_mean=float(y[-1]/th),
                           peak_cumulative_difference_from_endpoint_linear=float(max(abs(y))),
                           cumulative_half_gap=tip.tolist(),cumulative_full_gap=end.tolist(),
                           half_cumulative_crossings=fn(tip),full_cumulative_crossings=fn(end),
                           largest_energy_defect=max(abs(v['energy_defect']) for v in raw))
    return dict(scope='One physical aperture; test 50% of the archived hard moment replacement on two held unsplit grids.',
                exposure=x['history_exposure'],cases=c)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('depth1','depth4'):
        p.add_argument('--'+k,nargs=3,required=True,help='Native, half replacement, full replacement.')
    p.add_argument('--static',required=True)
    a=p.parse_args()
    def r(path):
        if path.endswith('.gz'):
            with gzip.open(path,'rt') as f:
                return json.load(f)
        return json.loads(Path(path).read_text())
    print(json.dumps(compare([r(f) for f in a.depth1],[r(f) for f in a.depth4],r(a.static)),indent=2,allow_nan=False))
