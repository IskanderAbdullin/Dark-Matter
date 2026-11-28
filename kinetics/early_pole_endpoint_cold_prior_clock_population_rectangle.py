"""Compare factors of earlier-clock-only cold-band response on exact same-chart panels."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_matched_clock_profile import fitted
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_spatial_population_exact_rectangle import decompose
from early_pole_endpoint_cold_spatial_population_factor_bounds import FIELDS,field
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def paired(pairs,cards,q,rates=None):
    defect=0.
    if not np.array_equal(cards[0][0].ds.ed,cards[1][0].ds.ed):
        raise ArithmeticError('Earlier timing cannot change spatial chart')
    local=[]
    for m in (0,1):
        books=[geometry(p[m]) for p in pairs]
        line={}
        for name,r,c in FIELDS:
            row=field(books,cards,q,r,c)
            b=np.asarray(row['native_and_crossed_bin_response_bounds'])
            if not np.array_equal(b[:,:,:,0],b[:,:,:,1]):
                raise ArithmeticError('An exact same-chart population substitution became interval-valued')
            v=decompose(b[:,:,:,0])
            if rates is not None and abs(np.expm1(v['native_log_ratio'])-rates[r][c])>3e-12:
                raise ArithmeticError('Native path does not close the paired receipt')
            line[name]=v
        local.append(line)
    for k,L in local[0].items():
        R=local[1][k]
        defect=max(defect,float(np.max(abs(np.asarray(L['response_rectangle'])/R['response_rectangle']-1.))),
                   float(np.max(abs(np.asarray(L['normalized_profile_rectangle'])-R['normalized_profile_rectangle']))))
    return local[1],defect


def run(pairs,states,archive,fixture):
    levels=(2,4,2,4)
    cards=[fitted(p,s,archive,n,2048) for p,s,n in zip(pairs,states,levels,strict=True)]
    if (fixture['state_hashes']!=[r[-1] for r in cards] or
        [fixture[k] for k in ('original_terminal_step_count','alternate_terminal_step_count')]!=[2048,2048]):
        raise ArithmeticError('Expected earlier-clock-only paired state evidence')
    q=np.asarray(fixture['common_q_edges'])
    if len(q)<6 or q[0]!=1e-6 or q[-1]!=1.5e-5:
        raise ArithmeticError('Need the exact reported common band')
    rows=[]
    defect=0.
    for i,n in enumerate((2,4)):
        I=[i+2,i] # Original earlier fine clock, then earlier altered coarse clock.
        v,err=paired([pairs[j] for j in I],[cards[j] for j in I],q,
                     fixture['earlier_coarse_over_original_band_minus_one'][i])
        rows.append(dict(scalar_subcells=n,observables=v))
        defect=max(defect,err)
    if defect>6e-5:
        raise ArithmeticError('Earlier-clock explicit population separation remains unresolved')
    return dict(x=.05,x_earlier_interval=[.02,.04],x_retained_interval=[.04,.05],
                original_and_alternate_earlier_steps=[2048,1024],matched_final_steps=2048,
                receipt_bound_state_hashes=fixture['state_hashes'],common_q_edges=q.tolist(),
                charts=rows,maximum_paired_relative_rate_or_absolute_profile_defect=defect,
                scope='Exact scalar-factor separation for one earlier-only clock perturbation; no decomposition of the held combined bare kernel.')

def run_one(n,pairs,states,archive,history,point,q):
    cards=[fitted(p,s,archive,n,2048) for p,s in zip(pairs,states,strict=True)]
    a,h=read(states[1])
    ref=cards[0][-1]
    if (n not in (2,4) or history['scalar_div']!=n or history['alternate_path']!=[1024,2048]
        or history['original_path']!=[2048,2048] or history['archive_prior_coarse_trace_defect']>5e-13
        or history['inputs_sha256'][2]!=ref or point['scalar_subcells']!=n
        or a['prior_clock_control']!='coarse_02_to_04'
        or point['original_and_alternate_state_sha256']!=[ref,h]
        or point['original_coarse_recovery_state_sha256']!=a['actual_prior_state_sha256']):
        raise ArithmeticError('The retained panels do not bind the earlier-only chronology')
    q=np.asarray(q,dtype=float)
    if len(q)!=18 or [q[0],q[-1]]!=[1e-6,1.5e-5]:
        raise ArithmeticError('Need the original 17-bin display')
    values,g=paired(pairs,cards,q)
    if g>6e-5:
        raise ArithmeticError('Same-chart four-corner response did not close')
    return dict(x=.05,scalar_subcells=n,common_q_edges=q.tolist(),earlier_interval=[.02,.04],
                original_and_alternate_earlier_steps=[2048,1024],matched_final_steps=2048,
                original_and_alternate_state_sha256=[z[-1] for z in cards],
                observables=values,maximum_paired_relative_rate_or_absolute_profile_defect=g,
                scope='Two earlier timing histories on one fixed scalar grid and identical final time-step count, not a physical production or time extrapolation.')



def cli_one(argv):
    p=argparse.ArgumentParser(description=run_one.__doc__)
    for k in ('original','alternate','states'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    for k in ('archive','history','point','display'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--level',type=int,choices=(2,4),required=True)
    a=p.parse_args(argv)
    js=lambda f: json.loads(f.read_text())
    point=next(row for row in js(a.point)['results'] if row['scalar_subcells']==a.level)
    return run_one(a.level,(a.original,a.alternate),a.states,a.archive,js(a.history),point,js(a.display)['common_q_edges'])


if __name__=='__main__':
    import sys
    if sys.argv[1:2]==['one']:
        print(json.dumps(cli_one(sys.argv[2:]),indent=2,allow_nan=False))
        sys.exit(0)
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s2-alt','s4-alt','s2-original','s4-original'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    p.add_argument('--states',type=Path,nargs=4,required=True)
    for k in ('archive','fixture'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run((a.s2_alt,a.s4_alt,a.s2_original,a.s4_original),a.states,a.archive,
                         json.loads(a.fixture.read_text())),indent=2,allow_nan=False))
