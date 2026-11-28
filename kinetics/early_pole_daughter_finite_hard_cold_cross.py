"""Measure a mixed finite discretization contrast rather than assuming separability."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(h1q8,h4q8,h1q32,h4q32,static):
    ref=h1q8[0]
    t=np.asarray(ref['history_exposure'])
    settings=('duration','projection','report_panels','rtol','base_atol','history_exposure')
    hist={}
    detail={}
    modes={}
    anchor=None
    for k,pair,depth,n in (('h1q8',h1q8,1,8),('h4q8',h4q8,4,8),
                            ('h1q32',h1q32,1,32),('h4q32',h4q32,4,32)):
        rows=[]
        energy=[]
        vals=[]
        for x,label in zip(pair,('native','local-daughter'),strict=True):
            if (any(x[f]!=ref[f] for f in settings) or x.get('aperture_subdivision',1)!=depth
                    or x.get('prepared_charged_swap','native')!=label
                    or x.get('selected_elastic_channel','both')!='both'
                    or not x.get('elastic_channels',True)):
                raise ValueError('Mismatched exposure, kinetics or intended intervention.')
            err=abs(x['energy_defect'])
            if not (err<3e-12 and abs(x['entropy_defect'])<3e-11 and abs(x['hard_number_defect'])<3e-12):
                raise ValueError('Original energy, entropy and hard-number gates failed.')
            energy.append(err)
            y=np.asarray(x['cold_count_history'])
            if y.shape!=t.shape or not np.all(np.isfinite(y)&(y>0.)):
                raise ValueError('Incomplete or inadmissible reported count.')
            if abs(y[0]-ref['cold_count_history'][0])>4e-18:
                raise ValueError('Different physically retained cold count.')
            rows.append(y-y[0])
            if n==32:
                q=np.geomspace(1e-6,1e-5,33)**3
                np.testing.assert_array_equal(np.cbrt((q[:-1]+q[1:])/2),x['cold_momentum_nodes_over_T'])
                m=np.asarray(x['cold_mode_count_history'])
                if m.shape!=(len(t),n) or not np.all(np.isfinite(m)&(m>0.)):
                    raise ValueError('Incomplete admissible 32-mode history.')
                if anchor is None:
                    anchor=m[0]
                np.testing.assert_array_equal(anchor,m[0])
                if np.max(abs(m.sum(axis=1)-y)/np.spacing(y))>4.:
                    raise ValueError('Mode history did not reproduce scalar sum.')
                vals.append(m-m[0])
            elif 'cold_mode_count_history' in x:
                if np.asarray(x['cold_mode_count_history']).shape!=(len(t),n):
                    raise ValueError('Wrong cold resolution for original mode mesh.')
        contrast=rows[1]-rows[0]
        hist[k]=contrast
        brackets=np.flatnonzero((contrast[:-1]<0.)&(contrast[1:]>=0.))
        detail[k]=dict(final_contrast=float(contrast[-1]/t[-1]),count_contrast=contrast.tolist(),
                       negative_to_nonnegative_brackets=[t[i:i+2].tolist() for i in brackets],
                       max_energy_error=max(energy))
        if n==32:
            modes[k]=vals[1]-vals[0]
        if depth==4:
            local=pair[1]['local_archive_support']
            v=static['cases']['refined'][str(n)]
            if len(local['indices'])!=8 or np.any(np.diff(local['indices'])!=1):
                raise ValueError('Original daughter interval did not retain its eight exact cells.')
            np.testing.assert_allclose(local['momentum_interval'],(1.4138888888888888,1.4666666666666666),
                                       rtol=0.,atol=2e-14)
            if not np.isclose(local['initialized_cold_rate'],v['local_rate'],rtol=0.,atol=2e-16):
                raise ValueError('Unexpected static local initial RHS.')
    mixed=(hist['h4q32']-hist['h4q8'])-(hist['h1q32']-hist['h1q8'])
    hard32=(modes['h4q32']-modes['h1q32'])[-1]
    return dict(scope='Daughter local-minus-native finite held contrast on a 2x2 discretization grid, not a continuum extrapolation.',
                exposure=t.tolist(),cases=detail,
                mixed_correction_to_final_mean=float(mixed[-1]/t[-1]),mixed_correction_by_exposure=mixed.tolist(),
                hard_refinement_effect_on_q32_modes=(hard32/t[-1]).tolist(),
                h4q32_positive_terminal_modes=int(np.count_nonzero(modes['h4q32'][-1]>0.)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('h1q8','h4q8','h1q32','h4q32'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local-daughter.')
    p.add_argument('--static',required=True)
    v=p.parse_args()
    r=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(*([r(f) for f in getattr(v,k)] for k in ('h1q8','h4q8','h1q32','h4q32')),
                             r(v.static)),indent=2,allow_nan=False))
