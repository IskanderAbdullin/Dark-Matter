"""Separate rule, represented scalar history and charge history at the actual warm endpoint."""
import argparse
import json
from itertools import product
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_point_factorial import factors
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_actual_source_reconstruction import model,sample
from early_pole_endpoint_warm_collision_volume_reference import WarmVolume
from early_pole_endpoint_warm_source_cells import local


def run(native, states, control, prefix, fine_reference):
    n,hs=read(native)
    hc,hashc=read(prefix)
    h,hashh=read(fine_reference)
    x=.05
    a,y=unpack(n,False,x)
    m,_,_,H=a.state(x,y)
    if (h['x'],h['native_sha256'],h['fixed_optics_mass'],h['fixed_clock_H'],h['volume_orders'])!=(x,hs,m,H,[96,128,160]):
        raise ArithmeticError('Changed the direct Volume160 witness')
    if hc['contrasts'][-1]['common_native_sha256']!=hs or hc['contrasts'][-1]['x']!=x or len(hc['mixed_trace'])!=5:
        raise ArithmeticError('Need the exact complete paired terminal histories')
    d0=load('fine')[0]
    store=[]
    for i,p in enumerate(states):
        c,f=read(p)
        u=control if i else None
        gauss(c,a,m,x,budget_control=u)
        tag,d,R=chart_for(c,budget_control=u)
        z=np.asarray(c['nodal_occupation'])
        expected=(hc['uniform_trace'][-1]['uncompressed_state_sha256'] if not i else
                  hc['mixed_trace'][-1]['state_sha256'])
        if (f!=expected or f!=h['paired_actual_results'][i]['input_sha256']
                or not np.array_equal(d.ds.ed,a.ds.ed) or
                any(not np.array_equal(getattr(d,k),getattr(a,k)) for k in ('wp','wd','qp','qd'))):
            raise ArithmeticError('Endpoint chart optics, geometries or hashes changed')
        store.append((d,R,z[d.ns:],model(d,z,R,m)))
    w=a.ds.ed[:-1]>=.01
    codes=['G7','mixed255','Volume160']
    outcomes={}
    rows={}
    for k in range(3):
        if k<2:
            d,R,*_=store[k]
        else:
            d=WarmVolume(d0,SPECS['fine'],160)
            R=d.repeats
        P=np.where(d.ds.p>=.01,d.ds.p**2/np.hypot(d.ds.p,m),0.)
        for i,j in product((0,1),repeat=2):
            z,err=sample(d,R,store[j][2],m,store[i][3])
            if k==2 and err>4e-12:
                raise ArithmeticError('Volume160 representation did not conserve the original per-cell moments')
            u,v=[local(d,z,P,x,m,H,a.ds.ed,o) for o in (4,8)]
            if np.any(v[~w]!=0.):
                raise ArithmeticError('Warm rate escaped its retained support')
            drift=float(np.max(abs(v-u)/np.maximum(v[:,0,None]+v[:,1,None],1e-100)))
            if drift>3e-6:
                raise ArithmeticError('Crossed chart radial stencil did not resolve')
            net=np.array([fsum(col) for col in v.T])
            if net.min()<=0.:
                raise ArithmeticError('No positive gross and net directional totals')
            outcomes[(k,i,j)]=net
            rows[(k,i,j)]=dict(rule=codes[k],scalar_moments=codes[i],charge_occupation=codes[j],
                             max_radial_gross_relative_defect=drift,sampled_moment_max_defect=err,
                             integrated_forward_inverse_net=net.tolist(),warm_cell_forward_inverse_net=v[w].tolist())
    for i,c in enumerate(h['paired_actual_results']):
        if abs(outcomes[(2,i,i)][2]/c['common_reference_scale']-1.)>3e-13:
            raise ArithmeticError('Crossed high-order diagonal failed independent direct reference')
        original=c['common_reference_scale']*(1+float(c['aggregate_errors'][0][0]))
        projected=c['common_reference_scale']*(1+float(c['aggregate_errors'][2][0]))
        if abs(outcomes[(i,i,i)][2]/projected-1.)>3e-13 or abs(original/projected-1.)>2e-9:
            raise ArithmeticError('A low-rule diagonal failed its independent actual/stored reference')
    low={k:v for k,v in outcomes.items() if k[0]<2}
    avg,mini=factors(low)
    n=outcomes[(2,0,0)]
    o=outcomes[(2,1,1)]
    s=outcomes[(2,1,0)]
    c=outcomes[(2,0,1)]
    scalar=(np.log(s/n)+np.log(o/c))/2
    charge=(np.log(c/n)+np.log(o/s))/2
    if np.max(abs(scalar+charge-np.log(o/n)))>2e-14:
        raise ArithmeticError('Higher-volume history effects did not add')
    return dict(x=x,common_native_state_sha256=hs,mixed_prefix_sha256=hashc,source_reference_sha256=hashh,
                column_names=['warm_forward','warm_inverse','warm_net'],
                low_cube_factor_order=['quadrature_node_rule','represented_scalar_history','charged_history'],
                low_cube_log_effects_mean=avg.tolist(),low_cube_log_effects_range=mini.tolist(),
                volume160_history_only_mean_log_effects=[scalar.tolist(),charge.tolist()],
                low_rule_mixed255_over_G7_diagonal=(outcomes[(1,1,1)]/outcomes[(0,0,0)]-1.).tolist(),
                high_rule_mixed_history_over_G7_history_diagonal=(o/n-1.).tolist(),
                cases=list(rows.values()),
                scope='Frozen matched terminal counterfactuals. Quadrature rule, represented scalar moments and stored charge populations are exchanged independently without re-evolving. No time-integrated physical causal claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','control','prefix','fine-reference'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--states',nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.states,a.control,a.prefix,a.fine_reference),indent=2,allow_nan=False))
