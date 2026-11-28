"""Independent all-pairs step integral for the discontinuous endpoint chart."""
import json
from math import fsum

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_decay_mode_clock_reference import G,KAPPA
from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def direct(d,ar,x,m):
    factor=np.longdouble
    x,m,mb=map(factor,(x,m,ar['mass_ratio']*x))
    q=np.asarray(d.ds.p,dtype=factor)
    Es=np.hypot(m,q)
    A=(x*x+m*m-mb*mb)/2
    D=np.sqrt((x*x-(m+mb)**2)*(x*x-(m-mb)**2))/2
    lo=(A*Es-D*q)/(m*m)
    hi=(A*Es+D*q)/(m*m)
    E_p=np.hypot(x,np.asarray(ar['endpoint_parent_momentum_edges'],dtype=factor))
    E_d=np.hypot(mb,np.asarray(ar['endpoint_fermion_momentum_edges'],dtype=factor))
    fs=np.asarray(ar['endpoint_scalar_occupations'])
    fd=np.asarray(ar['endpoint_fermion_occupations'])
    fp=np.asarray(ar['endpoint_parent_occupations'])
    pair_F=(1-fd[:,None])*fp[None,:]
    pair_B=fd[:,None]*(1-fp[None,:])
    scalar=[]
    for i,(l,h) in enumerate(zip(lo,hi,strict=True)):
        a=np.maximum(l,np.maximum(E_p[:-1][None,:],Es[i]+E_d[:-1,None]))
        b=np.minimum(h,np.minimum(E_p[1:][None,:],Es[i]+E_d[1:,None]))
        width=np.maximum(b-a,0.)
        phase=KAPPA*G*x*((x+mb)**2-m*m)/(Es[i]*q[i])*d.ds.w[i]
        scalar.append([float(phase*(1+fs[i])*np.sum(width*pair_F)),
                       float(phase*fs[i]*np.sum(width*pair_B))])
    return np.array(scalar)


def run():
    results=[]
    for g,s in SPECS.items():
        ar=load(g)[0]
        x=ar['x_end']
        d=CellChart(ar,s,1)
        m,_,_,H=d.state(x,d.initial)
        y=direct(d,ar,x,m)
        halves=[]
        for n in (1,2,4,8):
            d=CellChart(ar,s,n)
            ch,_=d.chart(x,m,4)
            _,f,b=event(ch,d.initial,d)
            z=np.c_[np.bincount(ch[3],f,minlength=d.ns),np.bincount(ch[3],b,minlength=d.ns)]
            er=np.abs(y-z)/(np.maximum(np.max(abs(y),axis=0),1e-40))
            defect=float(np.max(er))
            integrated=float(np.max(abs((z.sum(axis=0)-y.sum(axis=0))/y.sum(axis=0))))
            if max(defect,integrated)>4e-12:
                raise ArithmeticError(('Direct step aperture agrees neither per-mode nor globally',g,n,defect,integrated))
            halves.append(dict(subintervals_per_archived_bin=n,
                               per_mode_max_scaled_defect=defect,integrated_defect=integrated))
        results.append(dict(grid=g,actual_direct_forward_and_inverse=[fsum(y[:,0]),fsum(y[:,1])],
                            aperture_controls=halves))
    return results


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
