"""Identify which frozen native cells limit the warm collision-source rule."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_collision_volume_reference import WarmVolume,net
from early_pole_endpoint_warm_pressure_source_control import cell_weights


def local(d,y,P,x,m,H,edges,n):
    chart,defect=d.chart(x,m,n)
    c,A,B=event(chart,y,d)
    if abs(fsum(c*np.concatenate(d.energies(x,m))))>3e-13 or defect>2e-14:
        raise ArithmeticError('Lost collision energy conservation')
    if min(np.min(A),np.min(B))<0:
        raise ArithmeticError('Nonpositive gross collision weight')
    i=chart[3]
    group=np.searchsorted(edges,d.ds.p,side='right')-1
    if np.min(group)<0 or max(group)>=len(edges)-1:
        raise ArithmeticError('Some nodes lie outside their native support')
    coeff=np.array([np.bincount(group[i],weights=H*P[i]*r,minlength=len(edges)-1)
                    for r in (A,B,A-B)]).T
    if np.max(abs(coeff[:,0]-coeff[:,1]-coeff[:,2])/np.maximum(np.max(abs(coeff),axis=1),1e-100))>5e-13:
        raise ArithmeticError('Local forward/inverse/net ledger did not close')
    return coeff


def run(c):
    x=c['x_end']
    a,y=unpack(c,False,x)
    m,_,_,H=a.state(x,y)
    edges=a.ds.ed
    floor=float(min(edges[edges>=.01]))
    warm=edges[:-1]>=floor
    config=load('fine')[0],SPECS['fine']
    rows=[]
    def append(name,d,v,P):
        if (not np.array_equal(d.wp,a.wp) or not np.array_equal(d.wd,a.wd)
                or abs(fsum(d.ds.w*v[:d.ns])/fsum(a.ds.w*y[:a.ns])-1)>5e-14):
            raise ArithmeticError('Frozen control changed charge chart or scalar count')
        lo,hi=[local(d,v,P,x,m,H,edges,n) for n in (4,8)]
        err=np.max(abs(lo-hi)/np.maximum(hi[:,0,None]+hi[:,1,None],1e-100))
        if err>3e-6 or np.any(hi[~warm]!=0):
            raise ArithmeticError(('Unresolved warm-cell source or pressure outside warm support',err))
        rows.append((name,hi,err))
    append('W1',a,y,cell_weights(edges,m,floor))
    for n in (2,4):
        b=ScalarSubshellChart(*config,2,4,1e-5,warm_floor=.01,warm_div=n)
        append(f'W{n}',b,np.r_[y[np.searchsorted(edges,b.ds.p,side='right')-1],y[a.ns:]],
               cell_weights(b.ds.ed,m,floor))
    b=ScalarGauss2(*config,2,4,1e-5,gauss_floor=.01)
    p=b.ds.p
    append('Gauss2',b,b.transfer_step(y),np.where(p>=floor,p*p/np.hypot(p,m),0.))
    for n in (32,64):
        b=WarmVolume(*config,n)
        q=b.ds.p
        v=b.lift(y)
        P=np.where(q>=floor,q*q/np.hypot(q,m),0.)
        append(f'Volume{n}',b,v,P)
    ref=rows[-1][1]
    scale=fsum(ref[:,2])
    if scale<=0:
        raise ArithmeticError('Positive reference pressure source required')
    residual=rows[-2][1]-ref
    local_defect=np.max(abs(residual)/np.maximum(ref[:,0,None]+ref[:,1,None],1e-100),axis=1)
    worst=int(np.argmax(local_defect))
    envelope=fsum(abs(residual[:,2]))/scale
    if envelope>4e-5:
        raise ArithmeticError(('The total unsigned volume reference drift is unresolved',envelope))
    independent=net(b,v,[P],x,m,H)[0]
    if np.max(abs(np.array([fsum(t) for t in ref.T])-independent)/np.maximum(abs(independent),1e-100))>3e-13:
        raise ArithmeticError('Cell reduction fails independent total')
    cases=[]
    for name,vals,err in rows:
        exc=vals[:,2]-ref[:,2]
        i=np.argsort(-abs(exc))[0:8]
        cases.append(dict(rule=name,max_radial_mesh_defect=float(err),
                         summed_net_relative_error=fsum(exc)/scale,L1_residual_relative=fsum(abs(exc))/scale,
                         top_eight_cells=[dict(q=edges[j:j+2].tolist(),net_error_over_total_ref=exc[j]/scale,
                               ref_net_fraction=ref[j,2]/scale) for j in i],
                         warm_cell_forward_inverse_net=vals[warm].tolist()))
    return dict(x=x,mass_held=m,clock_held=H,warm_floor=floor,warm_native_edges=edges[np.searchsorted(edges,floor):].tolist(),
                ref_net_pressure=scale,mesh_reference=dict(global_L1_over_total=bool(envelope<=4e-5),
                    global_L1_fraction=envelope, strict_native_cell_gate_passed=bool(local_defect[worst]<=4e-5),
                    maximum_relative_to_gross_cell_scale=local_defect[worst],native_q_of_maximum=edges[worst:worst+2].tolist(),
                    corresponding_ref_net_fraction=ref[worst,2]/scale),cases=cases,
                scope='Frozen collision source at identical piecewise-constant occupations and fitted clock; not a time-integrated decomposition.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('S4',type=Path)
    print(json.dumps(run(read(p.parse_args().S4)),indent=2,allow_nan=False))
