"""Sharp affine source bounds on the actual finite nodes at fixed warm moments."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_source_cells import local


def hull(u,c,mean):
    """All feasible two-node mixtures give the extrema with two equalities."""
    u,c=map(np.asarray,(u,c))
    if (len(u)!=len(c) or len(u)<2 or np.min(np.diff(u))<=0.
            or not u[0]<=mean<=u[-1] or not np.all(np.isfinite(c))):
        raise ArithmeticError('Invalid local finite-node moment problem')
    sol=[]
    for i,j in zip(*np.where((u[:,None]<=mean)&(u[None,:]>=mean))):
        if i==j:
            sol.append((float(c[i]),int(i),int(j),0.))
        else:
            t=(mean-u[i])/(u[j]-u[i])
            sol.append((float((1-t)*c[i]+t*c[j]),int(i),int(j),float(t)))
    if not sol:
        raise ArithmeticError('No two-node convex witness')
    return min(sol),max(sol)


def take(d,y,repeats,x,m,H,n):
    geom,err=d.chart(x,m,n)
    if err>2e-14:
        raise ArithmeticError('Independent charge mapping failed')
    empty=np.r_[np.zeros(d.ns),y[d.ns:]]
    all_=np.r_[np.ones(d.ns),y[d.ns:]]
    _,A0,B0=event(geom,empty,d)
    _,A1,B1=event(geom,all_,d)
    if np.any(B0!=0.) or not np.array_equal(A1,2*A0):
        raise ArithmeticError('Scalar dependence is not the declared affine ledger')
    P=np.where(d.ds.p>=.01,d.ds.p**2/np.hypot(d.ds.p,m),0.)
    i=geom[3]
    unit=H*P[i]
    constant=np.bincount(i,weights=unit*A0,minlength=d.ns)
    slope=np.bincount(i,weights=unit*(A0-B1),minlength=d.ns)/d.ds.w
    pop=d.ds.w*y[:d.ns]
    predicted=fsum(constant+pop*slope)
    ref=fsum(local(d,y,P,x,m,H,d.ds.ed,n)[:,2])
    if ref<=0. or abs(predicted/ref-1.)>3e-13:
        raise ArithmeticError('The finite-node source did not reproduce independent assembly')
    rows=[]
    ofs=np.r_[0,np.cumsum(repeats)]
    for k,(a,b) in enumerate(zip(d.ds.ed[:-1],d.ds.ed[1:],strict=True)):
        if repeats[k]==1 or a<.01:
            continue
        j=slice(ofs[k],ofs[k+1])
        p=pop[j]; count=fsum(p)
        u=scaled_energy(d.ds.p[j],a,b,m)
        mean=fsum(p*u)/count
        lo,hi=hull(u,slope[j],mean)
        c=fsum(constant[j])
        vl,vh=[c+count*z[0] for z in (lo,hi)]
        obs=c+fsum(p*slope[j])
        if obs<vl-5e-14*abs(ref) or obs>vh+5e-14*abs(ref):
            raise ArithmeticError('Actual chart escaped independent finite-node hull')
        rows.append(dict(q=[a,b],warm_count=count,scaled_energy_mean=mean,
                         actual_net=obs,allowed_net=[vl,vh],witnesses=[list(lo[1:]),list(hi[1:])]))
    low,high=[fsum(z['allowed_net'][i] for z in rows) for i in (0,1)]
    if any(not np.isfinite(z) for z in (low,high)) or not 0.<low<=ref<=high:
        raise ArithmeticError('Unexpected finite warm hull positivity')
    return dict(net=ref,extrema=[low,high],local_rows=rows)


def run(native,g7,mixed,control,mixed_chain):
    n,sh=read(native)
    mc,hc=read(mixed_chain)
    r=mc['contrasts'][-1]
    x=.05
    if r['x']!=x or r['common_native_sha256']!=sh or len(mc['mixed_trace'])!=5:
        raise ArithmeticError('Need the actual verified terminal checkpoint')
    a,y=unpack(n,False,x)
    m,_,_,H=a.state(x,y)
    if r['fixed_native_mass']!=m:
        raise ArithmeticError('Changed held diagnostic mass')
    result=[]
    for P,sha,f in ((g7,mc['uniform_trace'][-1]['uncompressed_state_sha256'],None),
                    (mixed,mc['mixed_trace'][-1]['state_sha256'],control)):
        c,h=read(P)
        if h!=sha:
            raise ArithmeticError('Exact evolved node card changed')
        gauss(c,a,m,x,budget_control=f)
        tag,d,repeats=chart_for(c,budget_control=f)
        y=np.asarray(c['nodal_occupation'])
        lo,hi=[take(d,y,repeats,x,m,H,k) for k in (4,8)]
        mesh=max(abs(np.array([lo['net'],*lo['extrema']])/np.array([hi['net'],*hi['extrema']])-1.))
        if mesh>2e-6:
            raise ArithmeticError('Independent radial hull totals did not resolve')
        result.append(dict(chart=tag,input_sha256=h,relative_mesh_defect=mesh,**hi))
    a,b=result
    p,q=a['extrema'],b['extrema']
    rr=[q[0]/p[1]-1.,q[1]/p[0]-1.]
    actual=b['net']/a['net']-1.
    if not rr[0]-2e-14<=actual<=rr[1]+2e-14:
        raise ArithmeticError('Matched actual comparison outside the certified finite-node interval')
    return dict(x=x,native_sha256=sh,mixed_chain_sha256=hc,common_m=m,common_clock_H=H,
                paired_chart_results=result,mixed255_over_G7_relative_net=actual,
                allowed_relative_net_at_stored_node_moments=rr,
                scope='Node weights may vary freely at fixed count and mean energy in each warm cell, with the original nodes and charge state of each chart. These are conditional instantaneous affine sources, not continuum or dynamical extrema.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','g7','mixed','control','mixed-chain'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.g7,a.mixed,a.control,a.mixed_chain),indent=2,allow_nan=False))
