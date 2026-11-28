"""Disentangle discrete occupation from order error on actual terminal warm states."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import fit,scaled_energy
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_collision_volume_reference import WarmVolume,net
from early_pole_endpoint_warm_source_cells import local


def model(d,y,repeats,m):
    e=d.ds.ed
    ptr=np.r_[0,np.cumsum(repeats)]
    out=[]
    for i,(a,b) in enumerate(zip(e[:-1],e[1:],strict=True)):
        p=slice(ptr[i],ptr[i+1])
        N=d.ds.w[p]*y[p]
        s=fsum(N)
        if repeats[i]==1:
            out.append((None,float(y[p][0]),s,0.))
        else:
            v=scaled_energy(d.ds.p[p],a,b,m)
            mean=fsum(N*v)/s
            out.append((fit(a,b,m,s,mean),None,s,mean))
    return out


def sample(d,target_repeats,y,m,rows):
    e=d.ds.ed
    ptr=np.r_[0,np.cumsum(target_repeats)]
    occ=np.zeros(d.ns)
    maxerr=0.
    for i,(a,b) in enumerate(zip(e[:-1],e[1:],strict=True)):
        p=slice(ptr[i],ptr[i+1])
        card,c,n,u=rows[i]
        if card is None:
            occ[p]=c
            continue
        v=scaled_energy(d.ds.p[p],a,b,m)
        occ[p]=np.exp(card['log_occupation_at_lower_edge']+card['tilt']*v)
        w=d.ds.w[p]*occ[p]
        count=fsum(w)
        maxerr=max(maxerr,abs(count/n-1.),abs(fsum(w*v)/count-u))
    return np.r_[occ,y],maxerr


def run(native,g7,mixed,control,history,*,volumes=(32,64,96)):
    volumes=tuple(volumes)
    if volumes not in ((32,64,96),(96,128,160)):
        raise ArithmeticError('Use one independently declared three-volume ladder')
    base,sha=read(native)
    card,_=read(history)
    rr=card['contrasts'][-1]
    if (rr['x'],rr['common_native_sha256'])!=(.05,sha) or len(card['mixed_trace'])!=5:
        raise ArithmeticError('Use the two accepted actual terminal histories')
    a,y=unpack(base,False,.05)
    m,_,_,H=a.state(.05,y)
    if m!=rr['fixed_native_mass']:
        raise ArithmeticError('Reference clock and mass no longer match the fixed native snapshot')
    out=[]
    for p,truth,key,tag in ((g7,card['uniform_trace'][-1],'uncompressed_state_sha256','G7'),
                             (mixed,card['mixed_trace'][-1],'state_sha256','mixed255')):
        u,sha=read(p)
        f=control if tag=='mixed255' else None
        if sha!=truth[key]:
            raise ArithmeticError('Actual accepted input snapshot changed')
        gauss(u,a,m,.05,budget_control=f)
        label,d,repeats=chart_for(u,budget_control=f)
        y=np.asarray(u['nodal_occupation'])
        w=d.ds.ed[:-1]>=.01
        z=model(d,y,repeats,m)
        rec,err=sample(d,repeats,y[d.ns:],m,z)
        panels=[]
        mask=[]
        def keep(e,j,name,obs):
            if any(not np.array_equal(getattr(d,k),getattr(e,k)) for k in ('qp','qd','wp','wd')):
                raise ArithmeticError('The frozen control changed the charged nodes or projection')
            P=np.where(e.ds.p>=.01,e.ds.p**2/np.hypot(e.ds.p,m),0.)
            c,g=[local(e,j,P,.05,m,H,d.ds.ed,o) for o in (4,8)]
            drift=float(np.max(abs(c-g)/np.maximum(g[:,0,None]+g[:,1,None],1e-100)))
            v=net(e,j,[P],.05,m,H,8)[0]
            totals=np.array([fsum(k) for k in g.T])
            if drift>3e-6 or max(abs(totals/v-1.))>3e-13 or np.any(g[~w]!=0.):
                raise ArithmeticError('Warm-cell reference did not resolve its independent radial and total checks')
            mask.append(dict(name=name,radial_gross_defect=drift,warm_count_or_scaled_mean_error=obs))
            panels.append(g[w])
        keep(d,y,'actual_nodal',0.)
        keep(d,rec,'same_low_rule_sampled_reconstruction',err)
        for n in volumes:
            e=WarmVolume(load('fine')[0],SPECS['fine'],n)
            rec,err=sample(e,e.repeats,y[d.ns:],m,z)
            if err>4e-12:
                raise ArithmeticError('Fine reference failed the actual warm moments')
            keep(e,rec,f'volume_{n}_sampled_reconstruction',err)
        p,q,r,s,t=panels
        scale=fsum(t[:,2])
        if not scale>0.:
            raise ArithmeticError('Nonpositive warm source reference')
        diff=lambda u,v:(u[:,2]-v[:,2])/scale
        rdiff,sdiff=diff(r,s),diff(s,t)
        if max(fsum(abs(rdiff)),fsum(abs(sdiff)))>4e-5:
            raise ArithmeticError('Independent positive volumes did not resolve the warm net effect')
        closure,order=diff(p,q),diff(q,t)
        if np.max(abs(closure+order-diff(p,t)))>3e-15:
            raise ArithmeticError('Signed source components did not close independently')
        out.append(dict(chart=tag,input_sha256=sha,common_reference_scale=scale,checks=mask,
                        aggregate_errors=[[fsum(v),fsum(abs(v))] for v in (diff(p,t),closure,order,sdiff)],
                        all_cell_errors=[diff(p,t).tolist(),closure.tolist(),order.tolist(),sdiff.tolist()]))
    L,U=volumes[-2:]
    optional=dict(volume_orders=volumes) if volumes!=(32,64,96) else {}
    return dict(x=.05,native_sha256=rr['common_native_sha256'],fixed_optics_mass=m,fixed_clock_H=H,
                fields=[f'actual_minus_volume{U}','stored_vs_moment_fit_at_same_nodes',
                        f'same_fit_low_minus_volume{U}',f'volume{L}_minus_volume{U}'],**optional,
                scalar_cells=a.ds.ed[a.ds.ed>=.01].tolist(),paired_actual_results=out,
                scope='Terminal-state frozen source check, with the positive count/mean representative as a controlled reference, not a unique shape or integrated trajectory accuracy measure.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','g7','mixed','control','history'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--volumes',nargs=3,type=int,default=(32,64,96))
    a=p.parse_args()
    print(json.dumps(run(a.native,a.g7,a.mixed,a.control,a.history,volumes=a.volumes),indent=2,allow_nan=False))
