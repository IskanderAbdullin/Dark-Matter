"""Separate response update and time split effects for the aperture-local daughter."""
import argparse
import json
from pathlib import Path

import numpy as np
from early_pole_joined_temporal_refinement_reference import conservation


def compare(profiles, chart):
    conditions=('projection','grid','aperture_grid','hard_edge_sets','duration','pair','method',
                'strength','same_soft_strength','compton_strength','conversion_strength',
                'response_smear','elastic_energy_budget')
    a0,b0=profiles['active16']
    out={}
    for name,(a,b) in profiles.items():
        for ref,x,prep in ((a0,a,None),(b0,b,'local-daughter')):
            if any(x[k]!=ref[k] for k in conditions):
                raise ValueError(f'Mismatched split controls for {name}.')
            expected_steps=16 if name.endswith('16') else 64
            if x['steps']!=expected_steps or (x['screening_update']=='full')!=name.startswith('active'):
                raise ValueError(f'Unexpected split label, response or step count: {name}.')
            if prep and x.get('prepared_charged_swap')!=prep:
                raise ValueError(f'Unexpected split prepared input: {name}.')
            for k in ('cold_bins','hard_moments','instantaneous_source_by_q'):
                np.testing.assert_array_equal(x['snapshots'][0][k],ref['snapshots'][0][k])
        ca,cb=(np.asarray(x['state_history']) for x in (a,b))
        idx=np.arange(0,a['steps']+1,a['steps']//16)
        t=np.linspace(0,a['duration'],17)
        delta=((cb-cb[0])-(ca-ca[0]))[idx]
        out[name]=receipt(t,delta)
        out[name]['conservation']=[conservation(x) for x in (a,b)]
    a,b=chart
    if b.get('prepared_charged_swap')!='local-daughter' or any(
            a[k]!=b[k] for k in ('duration','projection','rtol','base_atol','history_exposure')):
        raise ValueError('Unexpected matched unsplit controls.')
    t=np.asarray(a['history_exposure'])
    np.testing.assert_array_equal(t,np.linspace(0,a0['duration'],17))
    ca,cb=(np.asarray(x['cold_count_history']) for x in (a,b))
    out['unsplit']=receipt(t,(cb-cb[0])-(ca-ca[0]))
    convergence={}
    if 'held64' in out:
        coarse=np.asarray(out['held16']['cumulative_effect'])
        fine=np.asarray(out['held64']['cumulative_effect'])
        exact=np.asarray(out['unsplit']['cumulative_effect'])
        # A diagnostic cancellation of a leading width-linear term, not a proven order.
        lead=(4*fine-coarse)/3
        convergence=dict(
            coarse_signed_mean_error=float((coarse[-1]-exact[-1])/t[-1]),
            fine_signed_mean_error=float((fine[-1]-exact[-1])/t[-1]),
            extrapolated_signed_mean_error=float((lead[-1]-exact[-1])/t[-1]),
            coarse_max_count_error=float(np.max(abs(coarse-exact))),
            fine_max_count_error=float(np.max(abs(fine-exact))),
            extrapolated_max_count_error=float(np.max(abs(lead-exact))),
            width_linear_mean_estimate=float(lead[-1]/t[-1]))
    return dict(cases=out, width_comparison=convergence,
                scope='Matched local daughter inputs: changes in the first-window integrated source.')


def receipt(t,d):
    if len(t)!=len(d):
        raise ValueError('Require matched history and exposure.')
    # Adjacent saved steps give a bracket; no interpolation claim.
    spans=[[float(t[i]),float(t[i+1])] for i in range(1,len(t)-1)
           if d[i]<0. and d[i+1]>=0.]
    return dict(final_mean_shift=float(d[-1]/t[-1]),
                cumulative_effect=d.tolist(),coarse_exposure=t.tolist(),
                running_mean_effect=(d[1:]/t[1:]).tolist(),
                negative_to_nonnegative_brackets=spans)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('active16','held16','held64','chart'):
        p.add_argument('--'+k,nargs=2,required=(k!='held64'),
                       help='Native and local-daughter output, in order.')
    q=p.parse_args()
    def read(paths):
        return [json.loads(Path(f).read_text()) for f in paths]
    seq={k:read(v) for k,v in vars(q).items() if k!='chart' and v}
    print(json.dumps(compare(seq,read(q.chart)),indent=2,allow_nan=False))
