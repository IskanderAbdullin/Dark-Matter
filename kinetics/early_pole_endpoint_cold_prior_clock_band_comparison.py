"""Compare true final bands while holding the final clock and varying the earlier clock."""
import argparse
import json
from pathlib import Path
from math import fsum

import numpy as np

from early_pole_endpoint_cold_matched_clock_profile import at_clock,common_partition,fitted
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def run(pairs,states,archive,history,point_receipt):
    if len(states)!=4 or len(pairs)!=4:
        raise ArithmeticError('Four actual terminal states are required')
    cards=[fitted(p,s,archive,level,2048) for p,s,level in zip(pairs,states,(2,4,2,4),strict=True)]
    for i,n in enumerate((2,4)):
        alt,h=read(states[i])
        ref=cards[i+2][-1]
        z=history[i]
        point=point_receipt['results'][i]
        if (z['scalar_div']!=n or z['alternate_path']!=[1024,2048] or z['original_path']!=[2048,2048]
            or z['archive_prior_coarse_trace_defect']>5e-13 or z['inputs_sha256'][2]!=ref
            or alt['prior_clock_control']!='coarse_02_to_04' or
            point['original_and_alternate_state_sha256']!=[ref,h]
            or point['original_coarse_recovery_state_sha256']!=alt['actual_prior_state_sha256']):
            raise ArithmeticError('Band does not bind the original and altered earlier clock records')
    out=[]
    for i in (0,1):
        b=[geometry(p[i]) for p in pairs]
        q=common_partition(b)
        profiles=[at_clock(z,w,q,independent_bracket=True) for z,w in ((b[:2],cards[:2]),(b[2:],cards[2:]))]
        whole=[np.array([[fsum(z[:,k,j]) for j in range(3)] for k in range(3)]) for _,z in b]
        shifts=(np.array(whole[:2])/whole[2:]-1).tolist()
        out.append((q.tolist(),profiles,shifts))
    (a,b,c),(d,e,f)=out
    err=float(np.max(abs(np.asarray(c)-f)))
    for u,v in zip(b,e,strict=True):
        for k in u:
            err=max(err,abs(u[k]['original_S2_S4_TV']-v[k]['original_S2_S4_TV']))
            err=max(err,abs(u[k]['count_and_mean_excess']['minimum']-v[k]['count_and_mean_excess']['minimum']))
    if a!=d or err>2e-9:
        raise ArithmeticError('Independent earlier-clock meshes do not reproduce the same conclusion')
    return dict(x=.05,k=.0015,common_q_edges=d,original_terminal_step_count=2048,alternate_terminal_step_count=2048,
                state_hashes=[x[-1] for x in cards],
                earlier_coarse_S2_S4_profiles=e[0],earlier_fine_S2_S4_profiles=e[1],
                earlier_coarse_over_original_band_minus_one=f,maximum_absolute_pair_defect=err,
                scope='Absolute and normalized finite response, altered x=.02–.04 timing only; not other intervals, continuum or physical error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s2_alt','s4_alt','s2_original','s4_original'):
        p.add_argument('--'+k.replace('_','-'),nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--history',nargs=2,type=Path,required=True)
    p.add_argument('--point-receipt',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run((a.s2_alt,a.s4_alt,a.s2_original,a.s4_original),a.states,a.archive,
                         [json.loads(p.read_text()) for p in a.history],json.loads(a.point_receipt.read_text())),indent=2,allow_nan=False))
