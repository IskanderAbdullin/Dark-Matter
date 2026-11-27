"""Compare exact local swaps on the same physical interval at two hard depths."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(original,refined,initial,global_continuous=None):
    collect={}
    whole={}
    ref=original[0]
    for profile,pair,depth in (('original',original,1),('refined',refined,4)):
        native,local=pair
        for tag,x in zip(('native','local-daughter'),pair,strict=True):
            controls=('duration','projection','report_panels','base_atol','rtol','history_exposure')
            if any(x[k]!=ref[k] for k in controls) or x.get('prepared_charged_swap','native')!=tag:
                raise ValueError('Mismatched chart, resolution or intervention.')
            if x.get('aperture_subdivision',1)!=depth or not x.get('elastic_channels',True) or (
                    x.get('selected_elastic_channel','both')!='both'):
                raise ValueError('Incorrect hard aperture or elastic terms.')
            if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                    and abs(x['hard_number_defect'])<3e-12):
                raise ValueError('Unaccepted original terminal gates.')
        memo=initial[profile]
        if len(memo['indices'])!=2*depth:
            raise ValueError('Wrong number of archived daughter cells.')
        if depth==4:
            meta=local['local_archive_support']
            if (meta['indices']!=memo['indices'] or meta['momentum_interval']!=memo['momentum_interval']
                    or meta['initialized_cold_rate']!=memo['local_initialized_cold_rate']):
                raise ValueError('The refined local chart does not match its static initialization.')
        a,b=(np.asarray(x['cold_count_history']) for x in pair)
        if (a.shape!=b.shape or len(a)!=ref['report_panels']+1 or a[0]!=b[0]
                or a[0]!=ref['cold_count_history'][0] or not np.all(np.isfinite((a,b)))):
            raise ValueError('Unexpected count history.')
        diff=b-a
        t=np.asarray(ref['history_exposure'])
        collect[profile]=dict(initial_rate_difference=memo['local_minus_native_initialized_rate'],
                              final_mean_local_shift=float(diff[-1]/ref['duration']),
                              cumulative_effect=diff.tolist(),
                              negative_to_nonnegative_saved_brackets=[
                                  t[i:i+2].tolist() for i in range(1,len(t)-1) if diff[i]<0<=diff[i+1]],
                              maximum_abs_energy_defect=max(abs(x['energy_defect']) for x in pair),
                              maximum_abs_entropy_defect=max(abs(x['entropy_defect']) for x in pair),
                              maximum_abs_hard_number_defect=max(abs(x['hard_number_defect']) for x in pair))
        if global_continuous is not None:
            v=global_continuous[profile]
            if (any(v[k]!=ref[k] for k in controls if k!='projection')
                    or v['projection']!='cell-conservative'
                    or v.get('prepared_charged_swap','native')!='native'
                    or v.get('aperture_subdivision',1)!=depth
                    or v.get('selected_elastic_channel','both')!='both'
                    or not v.get('elastic_channels',True)):
                raise ValueError('The whole archive control is different beyond its projection.')
            if not (abs(v['energy_defect'])<3e-12 and abs(v['entropy_defect'])<3e-11
                    and abs(v['hard_number_defect'])<3e-12):
                raise ValueError('Global continuous control did not pass original gates.')
            g=np.asarray(v['cold_count_history'])
            if g.shape!=a.shape or g[0]!=a[0] or not np.all(np.isfinite(g)):
                raise ValueError('The whole archive control changed the cold seed or history.')
            combined=g-a
            global_shift=float(combined[-1]/ref['duration'])
            whole[profile]=dict(global_final_mean_shift=global_shift,
                               local_over_global_finite_shift=float((b[-1]-a[-1])/combined[-1]) if combined[-1]!=0. else None,
                               whole_minus_local_terminal_mean=float((g[-1]-b[-1])/ref['duration']),
                               whole_minus_local_cumulative=(g-b).tolist())
    return dict(cases=collect,exposure=ref['history_exposure'],
                **({'whole_archive_context':whole} if whole else {}),
                finite_change=float(collect['refined']['final_mean_local_shift']-
                                    collect['original']['final_mean_local_shift']),
                scope='Only two projected geometry resolutions at a shared original momentum interval.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('original','refined'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then aperture-local.')
    p.add_argument('--global-continuous',nargs=2,help='Original then refined all-populations continuous controls.')
    p.add_argument('--static',required=True)
    a=p.parse_args()
    def read(f):
        return json.loads(Path(f).read_text())
    print(json.dumps(compare(list(map(read,a.original)),list(map(read,a.refined)),read(a.static),
                             (dict(zip(('original','refined'),map(read,a.global_continuous),strict=True))
                              if a.global_continuous else None)),
                     indent=2,allow_nan=False))
