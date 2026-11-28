"""Independently integrate the dominant warm error cell on exact charge-event cuts."""
import argparse
import hashlib
import os
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_actual_source_reconstruction import model
from early_pole_endpoint_warm_source_cells import local


def knots(d,x,m,a,b):
    mb=d.ratio*x
    pe,de=d.p_edges,d.d_edges
    P,F=np.hypot(pe,x),np.hypot(de,mb)
    A=(x*x+m*m-mb*mb)/2
    B=(x*x-m*m-mb*mb)/2
    D=np.sqrt((x*x-(m+mb)**2)*(x*x-(m-mb)**2))/2
    if min(A,B,D,m,mb)<=0. or a<=0. or b<=a:
        raise ArithmeticError('Invalid finite positive charge and scalar kinematics')
    # For a fixed charged endpoint, the scalar energies at the collinear
    # limits are (C E +/- D p)/m_charge^2. Rationalize the minus arm.
    u=(A*P+D*pe)/(x*x)
    l=(A*A+m*m*pe*pe)/(A*P+D*pe)
    v=(B*F+D*de)/(mb*mb)
    k=(B*B+m*m*de*de)/(B*F+D*de)
    e=np.r_[(P[:,None]-F).ravel(),u,l,v,k]
    e=e[(e>np.hypot(a,m))&(e<np.hypot(b,m))]
    q=np.sqrt((e-m)*(e+m))
    q=np.unique(np.r_[a,q[(q>a)&(q<b)],b])
    if len(q)<10 or np.min(np.diff(q))<=0.:
        raise ArithmeticError('Expected a nontrivial exact event topology')
    return q


def integral(c,edges,m,H,x,charge,fit,ns,nr,*,checkpoint_dir=None):
    u,w=roots_legendre(ns)
    a,b=edges[0],edges[-1]
    q=(edges[:-1,None]+np.diff(edges)[:,None]*(u+1)/2).ravel()
    v=(np.diff(edges)[:,None]*w/2).ravel()*q*q/(2*np.pi**2)
    if abs(fsum(v)/((b-a)*(a*a+a*b+b*b)/(6*np.pi**2))-1)>4e-14:
        raise ArithmeticError('Scalar volume factor changed across exact cuts')
    S=np.exp(fit['log_occupation_at_lower_edge']+fit['tilt']*scaled_energy(q,a,b,m))
    cache=None
    if checkpoint_dir is not None:
        key=hashlib.sha256()
        for u in (edges,np.asarray([m,H,x,ns,nr]),charge):
            key.update(np.ascontiguousarray(u,dtype='<f8').tobytes())
        key.update(json.dumps(fit,sort_keys=True,allow_nan=False).encode())
        key.update(bytes.fromhex('127a09ad1f2d8251f3a5b0fb8d3c7da9a76ae16ba572c9d83d779ed1882b7fbb'))  # unchanged integral kernel used for the persisted batches
        cache=checkpoint_dir/key.hexdigest()
        cache.mkdir(parents=True,exist_ok=True)
    rows=[]
    for i in range(0,len(q),128):
        dst=None if cache is None else cache/f'{i:06d}.json'
        if dst is not None and dst.exists():
            check=json.loads(dst.read_text())
            if check['start']!=i or check['stop']!=min(i+128,len(q)) or len(check['rate'])!=3 or not np.all(np.isfinite(check['rate'])):
                raise ArithmeticError('Persisted shell batch does not correspond to the exact panel selection')
            rows.append(check['rate'])
            continue
        p=q[i:i+128]
        d=ScalarSubshellChart(c,SPECS['fine'],2,4,1e-5)
        d.ds.ed=np.array([a,b])
        d.ds.p=p
        d.ds.w=v[i:i+128]
        d.ds.n=d.ns=len(p)
        d.iP=d.ns+d.nf
        d.n=d.iP+d.np
        d.w=np.r_[d.ds.w,d.wd,d.wp]
        occup=np.r_[S[i:i+128],charge]
        W=p*p/np.hypot(p,m)
        value=local(d,occup,W,x,m,H,d.ds.ed,nr)[0]
        if dst is not None:
            tmp=dst.with_suffix('.tmp')
            with tmp.open('w') as o:
                json.dump(dict(start=i,stop=min(i+128,len(q)),rate=value.tolist()),o,allow_nan=False)
                o.flush()
                os.fsync(o.fileno())
            tmp.replace(dst)
        rows.append(value)
    z=np.array([fsum(j) for j in np.asarray(rows).T])
    if z.min()<=0. or abs(z[0]-z[1]-z[2])>4e-13*z[0]:
        raise ArithmeticError('The batched scalar shell response does not close')
    return z,fsum(S*v)


def run(native,states,control,ref,cubes,*,checkpoint_dir=None):
    c,sha=read(native)
    r,hashr=read(ref)
    u,hashu=read(cubes)
    a,y=unpack(c,False,.05)
    m,_,_,H=a.state(.05,y)
    if sha!=r['native_sha256'] or u['source_reference_sha256']!=hashr or r['volume_orders']!=[96,128,160]:
        raise ArithmeticError('Do not silently select a different source snapshot or reference')
    labels=r['scalar_cells']
    ranks=[int(np.argmax(abs(np.asarray(d['all_cell_errors'][0])))) for d in r['paired_actual_results']]
    if len(set(ranks))!=1:
        raise ArithmeticError('Need the same dominant error cell in both actual stored rules')
    j=ranks[0]
    lo,hi=labels[j:j+2]
    edge=knots(a,.05,m,lo,hi)
    out=[]
    for i,s in enumerate(states):
        d,hashd=read(s)
        row=r['paired_actual_results'][i]
        f=control if i else None
        if hashd!=row['input_sha256']:
            raise ArithmeticError('Terminal actual state hash changed')
        gauss(d,a,m,.05,budget_control=f)
        tag,q,R=chart_for(d,budget_control=f)
        v=np.asarray(d['nodal_occupation'])
        index=np.flatnonzero(a.ds.ed[:-1]==lo)[0]
        z=model(q,v,R,m)[index]
        u1,u2=[integral(load('fine')[0],edge,m,H,.05,v[q.ns:],z[0],j,k,checkpoint_dir=checkpoint_dir)
               for j,k in ((4,8),(6,4))]
        if max(abs(u1[0]/u2[0]-1))>4e-6 or max(abs(j[1]/z[2]-1.) for j in (u1,u2))>3e-12:
            raise ArithmeticError('Independent exact-panel source rule or reconstructed count failed')
        tag=row['chart']
        former,=[j['warm_cell_forward_inverse_net'][j0] for j0 in (ranks[0],) for j in u['cases']
                 if (j['rule'],j['scalar_moments'],j['charge_occupation'])==('Volume160',tag,tag)]
        error=(np.array(former)-u2[0])/row['common_reference_scale']
        out.append(dict(chart=tag,original_positive_volume_160=former,direct_two_orders=[u1[0].tolist(),u2[0].tolist()],
                        measured_componentwise_order_difference=(u1[0]/u2[0]-1).tolist(),
                        volume160_minus_exact_shell_over_total_reference=error.tolist(),
                        actual_minus_exact_shell_net_over_total_reference=(row['all_cell_errors'][0][j]+error[2])))
    return dict(q_cell=[lo,hi],number_of_exact_scalar_panels=len(edge)-1,verified_source_cube_sha256=hashu,
                orders=[[4,8],[6,4]],shell_edges=edge.tolist(),charts=out,
                scope='One prespecified dominant terminal warm source cell. Exact topology cuts with two smooth-panel quadrature orders, not a whole-chart or physical error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','control','reference','cube'):
        p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--states',nargs=2,type=Path,required=True)
    p.add_argument('--checkpoint-dir',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.native,a.states,a.control,a.reference,a.cube,checkpoint_dir=a.checkpoint_dir),indent=2,allow_nan=False))
