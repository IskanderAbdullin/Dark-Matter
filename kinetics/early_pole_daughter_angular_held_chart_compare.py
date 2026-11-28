"""Conservation and like-for-like signed contrasts for an actual angular scalar run."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np


def compare(angular,control,expected):
    left,right=angular
    old,new=control
    ref=old
    t=np.asarray(ref['history_exposure'])
    settings=('projection','duration','report_panels','rtol','base_atol','history_exposure','aperture_subdivision')
    hist={}
    receipts={}
    for k,pair in (('ordinary',control),('angular',angular)):
        v=[]
        for arm,x in zip(('native','local-daughter'),pair,strict=True):
            if (any(x.get(a)!=ref.get(a) for a in settings) or x.get('prepared_charged_swap','native')!=arm
                    or not x.get('elastic_channels',True) or x.get('selected_elastic_channel','both')!='both'):
                raise ValueError('Different matched time chart or selected local intervention.')
            e,n,s=(abs(x[xkey]) for xkey in ('energy_defect','hard_number_defect','entropy_defect'))
            if not (e<3e-12 and n<3e-12 and s<3e-11):
                raise ValueError('An original solver gate did not pass.')
            a=np.asarray(x['cold_count_history'])
            if a.shape!=t.shape or not np.all(np.isfinite(a)&(a>0.)) or abs(a[0]-ref['cold_count_history'][0])>7e-18:
                raise ValueError('Inadmissible or physically distinct cold history.')
            if k=='angular':
                np.testing.assert_array_equal(x['cold_momentum_nodes_over_T'],expected['cold_momentum_nodes'])
                np.testing.assert_array_equal(x['initialized_scalar_volume_weights'],expected['scalar_volume_weights'])
                g=np.asarray(x['cold_mode_count_history'])
                if g.shape!=(len(a),expected['cold_nodes_count']) or not np.all(np.isfinite(g)&(g>0.)):
                    raise ValueError('Unexpected angular mode history.')
                np.testing.assert_array_equal(g[0],expected['seed'])
                # Match the storage layout used by the chart when it summed
                # modes before writing each row in transposed form.
                np.testing.assert_array_equal(g.T.copy(order='C').sum(axis=0),a)
                closure=float(np.max(abs(np.array([fsum(row) for row in g])-a)/np.spacing(a)))
                if not 0.<=closure<=8.:
                    raise ValueError('Unexpected compensated scalar closure defect.')
                receipts[arm+'_max_compensated_sum_difference_in_ulp']=closure
                v.append(g-g[0])
        hist[k]=(np.asarray(pair[1]['cold_count_history'])-pair[1]['cold_count_history'][0])-(np.asarray(pair[0]['cold_count_history'])-pair[0]['cold_count_history'][0])
        if k=='angular':
            h=v[1]-v[0]
            groups=np.asarray(expected['original_group'])
            grouped=np.stack([h[:,groups==i].sum(axis=1) for i in range(8)],axis=1)
            receipts['angular_grouped_signed_gains']=grouped.tolist()
            receipts['angular_terminal_positive_modes']=int(np.count_nonzero(h[-1]>0.))
            receipts['angular_terminal_negative_modes']=int(np.count_nonzero(h[-1]<0.))
    if left['angular_mesh_sha256']!=right['angular_mesh_sha256']:
        raise ValueError('Two angular meshes were not identical.')
    crossings={k:[t[i:i+2].tolist() for i in np.flatnonzero((x[:-1]<0)&(x[1:]>=0))]
               for k,x in hist.items()}
    return dict(scope='Only scalar sampling changes on a fixed held hard geometry. Not translatable to a physical relic.',
                angular_mesh_sha256=left['angular_mesh_sha256'],exposure=t.tolist(),
                negative_to_nonnegative_brackets=crossings,paired_means={k:float(x[-1]/t[-1]) for k,x in hist.items()},
                paired_count_histories={k:x.tolist() for k,x in hist.items()},**receipts)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--angular',nargs=2,required=True)
    p.add_argument('--ordinary',nargs=2,required=True)
    p.add_argument('--initializer',required=True)
    a=p.parse_args()
    def r(k):
        if k.endswith('.gz'):
            with gzip.open(k,'rt') as f:
                return json.load(f)
        return json.loads(Path(k).read_text())
    ar=[r(x) for x in a.angular]
    d=ar[0].get('aperture_subdivision',1)
    ex=next(p for p in r(a.initializer) if p['depth']==d)
    print(json.dumps(compare(ar,[r(x) for x in a.ordinary],ex),indent=2,allow_nan=False))
