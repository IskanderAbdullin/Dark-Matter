"""Integrated fixed-pole cold response on an opt-in continuum endpoint reconstruction."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_fixed_pole_scalar_q_band_reference import q_knots
from finite_mass_qed_timelike_pole_angular_reference import angle


def run(card,mesh,output,topology_only=False):
    g=card['grid']
    bundle=load(g)[0]
    tags={}
    if 'scalar_div' in card:
        from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
        d=ScalarSubshellChart(bundle,SPECS[g],card['charge_div'],card['scalar_div'],card['cut'])
        bundle=dict(bundle,scalar_momentum_edges=d.ds.ed)
        serial=json.dumps(card,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        tags=dict(scalar_div=card['scalar_div'],charge_div=card['charge_div'],scalar_cut=card['cut'],
                  state_sha256=hashlib.sha256(serial).hexdigest())
    else:
        d=CellChart(bundle,SPECS[g],card['subcells_per_archived_cell'],reference_x=card['reference_x'])
    y=np.asarray(card['endpoint_occupation'])
    if card['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid endpoint occupation card')
    x=card['x_end']
    m=d.state(x,y)[0]
    C=[]; fit=[]; callbacks=[]
    for mass,edges,n in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,err,s,alpha,beta=continuum_moments(edges,mass,n,36)
        fit.append((p,w,mass))
        callbacks.append(model_callback(mass,edges,alpha,beta))
        C.append(piecewise_logit(mass,edges,alpha,beta,n,1))
    k=.0015
    _,omega,res=root(lambda b:sum(charge_field(*a,b,36) for a in fit),k)
    span=(1e-6,1.5e-5)
    if d.ds.ed[0]>span[0] or d.ds.ed[-1]<span[1]:
        raise ArithmeticError('Missing retained scalar occupation support')
    u,v=[q_knots(bundle,C,k,omega,m,span,s) for s in (49,89)]
    if u.shape!=v.shape or max(abs(u-v))>9e-16:
        raise ArithmeticError('Shifted cold topology not resolved')
    meshes=((2,3,6,12),(4,6,10,18),(8,10,16,24))
    order,angular,inn,az=meshes[mesh]
    xs,ws=roots_legendre(order)
    expected={'grid':g,'time_x':x,'scalar_mass':m,'omega':omega,'residue':res,'mesh':list(meshes[mesh]),'total':len(v)-1}
    expected.update(tags)
    if topology_only:
        return dict(expected,edges=v.tolist())
    done=set()
    if output.exists():
        for line in output.read_text().splitlines():
            r=json.loads(line)
            j=r['panel']
            if j in done or r['bounds']!=v[j:j+2].tolist() or any(r[k]!=a for k,a in expected.items()):
                raise ArithmeticError('Changed endpoint or mesh in resume')
            done.add(j)
    with output.open('a',encoding='utf-8',buffering=1) as f:
        for j,(a,b) in enumerate(zip(v[:-1],v[1:],strict=True)):
            if j in done:
                continue
            coeff=np.zeros((3,3))
            for q,z in zip((a+b)/2+(b-a)*xs/2,(b-a)*ws/2,strict=True):
                E=np.hypot(q,m)
                n=y[np.searchsorted(d.ds.ed,q,side='right')-1]
                A,B,D=angle(C,k,omega,res,E,'T',angular,(inn,az),q,occupation_callbacks=callbacks)[0]
                data=z*q*q/E*np.array([(1+n)*A,n*B,n*D+A])
                coeff[0]+=data
                coeff[1]+=E*data
                coeff[2]+=(q*q/(E+m))*data
            print(json.dumps(dict(expected,panel=j,bounds=[a,b],coeff=coeff.tolist()),allow_nan=False),
                  file=f,flush=True)
            os.fsync(f.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('mesh',type=int,choices=(0,1,2))
    p.add_argument('output',type=Path,nargs='?')
    p.add_argument('--topology-only',action='store_true')
    o=p.parse_args()
    if o.output is None and not o.topology_only:
        p.error('output is required except in topology-only mode')
    result=run(json.loads(o.state.read_text()),o.mesh,o.output,o.topology_only)
    if result is not None:
        print(json.dumps(result,indent=2,allow_nan=False))
