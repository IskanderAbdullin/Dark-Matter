"""Check paired per-mode cold contrasts after tighter integration."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(base,tight):
    ref=base[0]
    fields=('duration','projection','report_panels','history_exposure','cold_momentum_nodes_over_T')
    rows={}
    closure={}
    arms_at_tolerance={}
    seed=np.asarray(ref['cold_mode_count_history'][0])
    if len(seed)!=32 or not np.all(np.isfinite(seed)&(seed>0.)):
        raise ValueError('Require all positive initial archived counts.')
    u=np.geomspace(1e-6,1e-5,33)**3
    np.testing.assert_array_equal(np.cbrt((u[:-1]+u[1:])/2),
                                  ref['cold_momentum_nodes_over_T'])
    for tag,arms,rtol,atol in (('nominal',base,3e-9,2e-14),('tight',tight,3e-10,2e-15)):
        data=[]
        closure[tag]={}
        for label,x in zip(('native','local-daughter'),arms,strict=True):
            if any(x[k]!=ref[k] for k in fields) or x['rtol']!=rtol or x['base_atol']!=atol:
                raise ValueError('Mismatched chart data beyond the tolerance.')
            if (x.get('prepared_charged_swap','native')!=label or x.get('aperture_subdivision',1)!=1
                    or not x.get('elastic_channels',True) or x.get('selected_elastic_channel','both')!='both'):
                raise ValueError('Wrong local preparation or held kernel.')
            if not (abs(x['energy_defect'])<3e-12 and abs(x['entropy_defect'])<3e-11
                    and abs(x['hard_number_defect'])<3e-12):
                raise ValueError('Original gates not accepted.')
            v=np.asarray(x['cold_mode_count_history'])
            if v.shape!=(ref['report_panels']+1,32):
                raise ValueError('Expected all 32 mode histories.')
            if not np.all(np.isfinite(v)&(v>0.)):
                raise ValueError('A recorded cold count is not admissible.')
            totals=np.asarray(x['cold_count_history'])
            consistency=np.abs(v.sum(axis=1)-totals)/np.spacing(totals)
            if not np.all(np.isfinite(consistency)) or np.max(consistency)>4.:
                raise ValueError('Stored per-mode histories do not reproduce the scalar count.')
            closure[tag][label]=float(np.max(consistency))
            np.testing.assert_array_equal(v[0],seed)
            data.append(v)
        rows[tag]=data[1]-data[0]
        arms_at_tolerance[tag]=data
    a,b=(d[-1] for d in rows.values())
    abs_change=abs(b-a)
    ratio=np.divide(abs_change,abs(b),out=np.full_like(b,np.inf),where=b!=0.)
    t=ref['duration']
    signed_resolution=b/np.spacing(seed)
    if a.sum()==0.:
        raise ValueError('Unresolved nominal total paired increment.')
    rel_gap=(b.sum()-a.sum())/a.sum()
    sample=np.asarray(ref['history_exposure'])
    crossing_brackets={}
    for tag,contrast in rows.items():
        crossing_brackets[tag]=[[sample[i:i+2].tolist() for i in
                                 np.flatnonzero((one[:-1]<0.)&(one[1:]>=0.))]
                                for one in contrast.T]
    a0,b0=(crossing_brackets[tag] for tag in ('nominal','tight'))
    crossing_changes=[k for k,(x,y) in enumerate(zip(a0,b0,strict=True)) if x!=y]
    arm_gaps={}
    for j,label in enumerate(('native','local-daughter')):
        v=arms_at_tolerance['tight'][j]-arms_at_tolerance['nominal'][j]
        arm_gaps[label]=dict(nominal_rhs_calls=base[j]['calls'],tight_rhs_calls=tight[j]['calls'],
                             peak_abs_mode_count_change=float(np.max(np.abs(v))),
                             max_terminal_mode_count_change_in_initial_spacings=float(
                                 np.max(np.abs(v[-1])/np.spacing(seed))))
    return dict(terminal_means_from_modes={k:d[-1].sum()/t for k,d in rows.items()},
                terminal_paired_mean_by_mode={k:(d[-1]/t).tolist() for k,d in rows.items()},
                smallest_positive_terminal_count=float(min(b[b>0.])) if np.any(b>0.) else None,
                tight_terminal_count_in_initial_mode_spacings=signed_resolution.tolist(),
                minimum_tight_pair_over_initial_mode_spacing=float(np.min(signed_resolution)),
                maximum_sum_to_scalar_spacing_defect=closure,
                individual_tolerance_changes=arm_gaps,
                relative_change_in_paired_final_mean=float(rel_gap),
                saved_cumulative_crossing_brackets=crossing_brackets,
                modes_with_changed_cumulative_crossing_brackets=crossing_changes,
                modes_changing_or_losing_sign=np.flatnonzero(np.sign(a)!=np.sign(b)).tolist(),
                negative_or_zero_fine_mode_indices=np.flatnonzero(b<=0.).tolist(),
                maximum_absolute_change_in_paired_count_history=float(np.max(abs(rows['tight']-rows['nominal']))),
                maximum_terminal_change_over_tight_magnitude=float(np.max(ratio[np.isfinite(ratio)]))
                  if np.any(np.isfinite(ratio)) else None,
                scope='An empirical two-tolerance test on the same original-hard, 32-cold discretization.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('base','tight'):
        p.add_argument('--'+k,nargs=2,required=True,help='Native then local.')
    a=p.parse_args()
    print(json.dumps(compare(*([json.loads(Path(f).read_text()) for f in paths]
                                 for paths in vars(a).values())),indent=2,allow_nan=False))
