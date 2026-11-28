"""Conservative direction-wise capacity limiter on the frozen charged chart."""
import json
from math import fsum

import numpy as np

from early_pole_double_charge_stencil_pilot import chart,event
from early_pole_edge_padded_charged_projection_pilot import charge_mesh
from early_pole_endpoint_stencil_remap_pilot import LABELS,SPECS,archives,load,parts,remap,MassiveFermion
from early_pole_live_parent_energy_reference import state_moments


def step(ch,d,f,dt,*,rate_state=None):
    _,fo,inv=event(ch,f if rate_state is None else rate_state,d)
    extent=dt*(fo-inv)
    _,_,_,i,j,k,a,b,_=ch
    ix=(i,d.ns+j-1,d.ns+j,d.iP+k-1,d.iP+k)
    dep=(np.ones_like(a),1-a,a,-(1-b),-b)
    flows=[x*extent for x in dep]
    outgoing=sum(np.bincount(q,np.maximum(-s,0.),minlength=d.n) for q,s in zip(ix,flows,strict=True))
    incoming=sum(np.bincount(q,np.maximum(s,0.),minlength=d.n) for q,s in zip(ix,flows,strict=True))
    stored=f*d.w
    holes=np.r_[np.full(d.ns,np.inf),d.w[d.ns:]*(1-f[d.ns:])]
    slack=1-32*np.finfo(float).eps
    em=np.divide(slack*stored,outgoing,out=np.ones_like(f),where=outgoing>0)
    im=np.divide(slack*holes,incoming,out=np.ones_like(f),where=incoming>0)
    th=np.minimum.reduce([np.where(s<0,em[q],im[q]) for q,s in zip(ix,flows,strict=True)])
    th=np.clip(th,0.,1.)
    next_delta=sum(np.bincount(q,s*th,minlength=d.n) for q,s in zip(ix,flows,strict=True))
    new=f+next_delta/d.w
    if np.min(new)<-1e-16 or np.max(new[d.ns:])>1+1e-14:
        raise ArithmeticError('A limited event step left the admissible occupations')
    return new,float(np.sum(abs(extent)*(1-th))/max(np.sum(abs(extent)),1e-100))


def ssp(ch,d,f,dt):
    a,c=step(ch,d,f,dt)
    b,e=step(ch,d,a,dt)
    return f+.5*(b-f),max(c,e)


def run():
    rows=[]
    for g,shape in SPECS.items():
        a=load(g)[0]
        x=a['x_end']
        d=MassiveFermion(*shape,a['mass_ratio'],a['elastic_rate'],True)
        charge_mesh(d,1/4096,1/4096)
        y=np.concatenate([remap(o,a['endpoint_'+k+'_occupations'],n)[0]
                          for o,k,n in zip(archives(a),LABELS,parts(d),strict=True)])
        E,ra,rb,m=state_moments(d,x,y)
        H=np.sqrt(d.radiation/(d.radiation+ra+rb))
        ch,_=chart(d,x,m,4)
        for stepper in (step,ssp):
            splits=[]
            for n in (32,64,128,256,512):
                f=y.copy()
                cap=[]
                for _ in range(n):
                    f,c=stepper(ch,d,f,H*1e-6/n)
                    cap.append(c)
                w=(f-y)*d.w
                dn=[float(fsum(z)) for z in np.split(w,(d.ns,d.iP))]
                err=abs(float(fsum(w*E)))
                if min(f)<0 or max(f[d.ns:])>1 or max(abs(dn[0]-dn[1]),abs(dn[0]+dn[2]),err)>5e-14:
                    raise ArithmeticError('Conservative event capacities lost occupation, count or fixed energy')
                splits.append(dict(steps=n,delta_species=dn,energy_residual=float(fsum(w*E)),
                                   min_occupation=float(min(f)),max_fermion=float(max(f[d.ns:])),
                                   max_limited_absolute_extent_fraction=max(cap)))
            n=np.array([x['delta_species'][0] for x in splits])
            d1,d2=np.diff(n[-3:])
            if abs(d2/n[-1])>1.6e-4 or not 0<abs(d2)<abs(d1):
                raise ArithmeticError(('Frozen time refinement was not resolved',g,stepper.__name__,n))
            trial=dict(grid=g,scheme=stepper.__name__,scope='Frozen chart advance, no mass or expansion feedback',
                       last_scalar_number_step_relative=abs(d2/n[-1]),last_delta_ratio=float(d1/d2),
                       trials=splits)
            rows.append(trial)
    return rows


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
