"""Match cold profile resolution across two terminal clock refinements."""
import argparse
import json
from functools import reduce
from math import fsum
from pathlib import Path
import numpy as np

from early_pole_endpoint_cold_normalized_profile_envelope import build,geometry,compare
from early_pole_endpoint_cold_profile_exact_finite_distance import solve
from early_pole_endpoint_cold_scalar_band_ladder import verify_pair
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def fitted(paths,state,archive,level,steps):
    d=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,level,1e-5)
    c,sha=read(state)
    r=verify_pair(paths,state,archive,level,steps=steps)
    if sha!=r['state_sha256_uncompressed']:
        raise ArithmeticError('Unverified endpoint state')
    v=np.asarray(c['endpoint_occupation'])
    if v.shape!=(d.n,):
        raise ArithmeticError('Wrong chart for the stored endpoint')
    return d,v[:d.ns],r['diagnostic']['full_pole_and_snapshot']['scalar_mass'],sha


def common_partition(books):
    q=reduce(np.intersect1d,(z[0] for z in books))
    if len(q)<6 or [q[0],q[-1]]!=[1e-6,1.5e-5]:
        raise ArithmeticError('Insufficient exactly shared retained q support')
    return q


def at_clock(books,states,bins):
    r={}
    for obs in ('forward_number','inverse_number','forward_excess'):
        a,b=[build(p,s,bins,raw=True,observable=obs) for p,s in zip(books,states,strict=True)]
        m,n=[build(p,s,bins,observable=obs) for p,s in zip(books,states,strict=True)]
        orig=fsum(abs(a['base']-b['base']))/2
        u=compare(m,n,bins)
        row={'original_S2_S4_TV':orig,'S4_over_S2_absolute_band_response':b['total']/a['total']}
        for k in ('count_only','count_and_mean_excess'):
            s=solve(a,b,k)
            val=s['minimum_finite_common_bin_TV']
            lo=u[k]['robust_total_variation_lower_bound']
            if not lo-2e-9<=val<=orig+2e-9:
                raise ArithmeticError('Finite simplex result conflicts with known bounds')
            row[k]=dict(minimum=val,projected_profile_TV=s['direct_witness_TV'],regrouped_lower_bound=lo)
        r[obs]=row
    return r


def run(pairs,states,archive,clock_receipt):
    if len(states)!=4 or len(pairs)!=4:
        raise ValueError('Need same restart coarse and fine S2/S4 for four snapshots')
    c=[fitted(p,s,archive,n,t) for p,s,n,t in zip(pairs,states,(2,4,2,4),(1024,1024,2048,2048),strict=True)]
    # The integrated program orders its independent clock records by level,
    # rather than by matched clock. Verify exact state bytes before regrouping.
    if clock_receipt['charts']!=[2,2,4,4] or clock_receipt['step_counts']!=[1024,2048,1024,2048]:
        raise ArithmeticError('The integrated counterpart uses a different chart/clock schedule')
    ordered=[x['state_sha256_uncompressed'] for x in clock_receipt['full_state_panel_receipts']]
    if len(ordered)!=4 or [ordered[i] for i in (0,2,1,3)]!=[x[-1] for x in c]:
        raise ArithmeticError('The profile and integral artifacts use different state bytes')
    meshes=[]
    for i in (0,1):
        b=[geometry(p[i]) for p in pairs]
        q=common_partition(b)
        meshes.append((q.tolist(),[at_clock(b[:2],c[:2],q),at_clock(b[2:],c[2:],q)]))
    (lo,u),(hi,v)=meshes
    defects=[]
    for x,y in zip(u,v,strict=True):
        for k,a in x.items():
            b=y[k]
            defects.extend([abs(a['original_S2_S4_TV']-b['original_S2_S4_TV']),
                            abs(a['S4_over_S2_absolute_band_response']-b['S4_over_S2_absolute_band_response'])])
            defects.extend(abs(a[j]['minimum']-b[j]['minimum']) for j in ('count_only','count_and_mean_excess'))
    err=max(defects)
    if hi!=lo or err>2e-9:
        raise ArithmeticError(('Different comparison support or unresolved mesh',err))
    first,second=v
    shifts={}
    for k,a in first.items():
        b=second[k]
        shifts[k]={j:b[j]-a[j] for j in ('original_S2_S4_TV','S4_over_S2_absolute_band_response')}
        shifts[k].update({j:b[j]['minimum']-a[j]['minimum'] for j in ('count_only','count_and_mean_excess')})
    return dict(x=.05,common_q_edges=hi,ordered_state_sha256_uncompressed=[x[-1] for x in c],
                coarse=first,fine=second,fine_minus_coarse=shifts,maximum_pair_mesh_absolute_defect=err,
                scope='Discrete absolute and normalized cold responses at matched clocks, held optics and exact four-way q partition. No earlier-history or continuum guarantee.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('s2_coarse','s4_coarse','s2_fine','s4_fine'):
        p.add_argument('--'+k.replace('_','-'),nargs=2,type=Path,required=True)
    p.add_argument('--states',nargs=4,type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--integrated-clock-receipt',type=Path,required=True)
    a=p.parse_args()
    r=json.loads(a.integrated_clock_receipt.read_text())
    print(json.dumps(run((a.s2_coarse,a.s4_coarse,a.s2_fine,a.s4_fine),a.states,a.archive,r),indent=2,allow_nan=False))
