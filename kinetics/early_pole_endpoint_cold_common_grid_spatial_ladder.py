"""Compare adjacent scalar levels on one exact common momentum partition."""
import argparse
import json
from functools import reduce
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import build,fit_state,geometry
from early_pole_endpoint_cold_profile_exact_finite_distance import solve


def compare(pairs,states,levels,archive,*,independent_bracket=False):
    if levels not in ((2,4,8),(2,4,8,16)):
        raise ArithmeticError('Use the consecutive completed spatial ladder')
    S=[fit_state(p,archive,n,pair) for p,n,pair in zip(states,levels,pairs,strict=True)]
    out=[]
    for t in (0,1):
        G=[geometry(pair[t]) for pair in pairs]
        q=reduce(np.intersect1d,(g[0] for g in G))
        if len(q)<6 or q[0]!=1e-6 or q[-1]!=1.5e-5:
            raise ArithmeticError('The common comparison lost its fixed cold band')
        rows=[]
        for i in range(len(levels)-1):
            r={}
            for obs in ('forward_number','inverse_number','forward_excess'):
                a,b=[build(g,p,q,raw=True,observable=obs) for g,p in zip(G[i:i+2],S[i:i+2],strict=True)]
                signed=b['base']-a['base']
                j=int(np.argmax(abs(signed)))
                o=fsum(abs(signed))/2
                m={}
                for k in ('count_only','count_and_mean_excess'):
                    v=solve(a,b,k,independent_bracket=independent_bracket)
                    if not -2e-9<=v['minimum_finite_common_bin_TV']<=o+2e-9:
                        raise ArithmeticError('Known original profile is not a feasible upper witness')
                    m[k]=dict(minimum=v['minimum_finite_common_bin_TV'],direct_witness=v['direct_witness_TV'])
                    if independent_bracket:
                        cert=v['independent_profile_bracket']
                        m[k].update({key:cert[key] for key in ('support_lower_bound','projected_profile_upper_bound','verified_TV_bracket_width')})
                r[obs]=dict(finer_over_coarser_band_response=b['total']/a['total'],original_TV=o,largest_signed_bin_edges=q[j:j+2].tolist(),signed_profile_changes=signed.tolist(),**m)
            rows.append(dict(levels=list(levels[i:i+2]),observables=r))
        out.append((q.tolist(),rows))
    (u,x),(v,y)=out
    if u!=v:
        raise ArithmeticError('Independent q meshes disagree on the common partition')
    defect=0.
    for a,b in zip(x,y,strict=True):
        for k,t in a['observables'].items():
            w=b['observables'][k]
            d=[abs(t[j]-w[j]) for j in ('original_TV','finer_over_coarser_band_response')]
            d.extend(abs(t[j]['minimum']-w[j]['minimum']) for j in ('count_only','count_and_mean_excess'))
            defect=max(defect,*d)
    if defect>2e-9:
        raise ArithmeticError('Mesh comparison did not close')
    return dict(x=.05,k=.0015,charge_subcells=2,levels=levels,states_sha256_uncompressed=[z[-1] for z in S],
                q_edges=v,adjacent=y,maximum_paired_mesh_difference=defect,
                scope='Fixed held optical kernels and finite positive panel reweights at each native chart own moments; common coarse display only, not continuum or asymptotic convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in (2,4,8):
        p.add_argument(f'--s{k}',nargs=2,type=Path,required=True)
    p.add_argument('--s16',nargs=2,type=Path)
    p.add_argument('--states',nargs='+',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--independent-bracket',action='store_true')
    a=p.parse_args()
    q=(a.s2,a.s4,a.s8)+((a.s16,) if a.s16 else ())
    print(json.dumps(compare(q,a.states,(2,4,8,16)[:len(q)],a.archive,independent_bracket=a.independent_bracket),indent=2,allow_nan=False))
