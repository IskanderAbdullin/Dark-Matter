"""Compare several frozen warm scalar collision quadratures with a higher-order reference."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_double_charge_stencil_pilot import event
from early_pole_endpoint_scalar_gauss2_live_fixture import ScalarGauss2
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from early_pole_endpoint_warm_pressure_source_control import cell_weights


class WarmVolume(ScalarSubshellChart):
    def __init__(self,a,shape,o):
        super().__init__(a,shape,2,4,1e-5)
        edge=self.ds.ed
        lo,hi=edge[:-1],edge[1:]
        self.nold=self.ns
        self.take=lo>=.01
        t,w=roots_legendre(o)
        h=hi-lo
        q=lo[:,None]+h[:,None]*(t+1)/2
        weight=w*q*q*h[:,None]/(4*np.pi**2)
        weight*=self.ds.w[:,None]/weight.sum(axis=1)[:,None]
        if not np.allclose(weight.sum(axis=1),self.ds.w,rtol=4e-15,atol=0):
            raise ArithmeticError('Warm reference quadrature lost the scalar shell volume')
        self.repeats=np.where(self.take,o,1)
        sel=np.arange(o)[None,:]<self.repeats[:,None]
        q[~self.take,0]=self.ds.p[~self.take]
        weight[~self.take,0]=self.ds.w[~self.take]
        self.ds.p=q[sel]
        self.ds.w=weight[sel]
        self.ds.n=len(self.ds.w)
        self.ns=self.ds.n
        self.iP=self.ns+self.nf
        self.n=self.iP+self.np
        self.w=np.r_[self.ds.w,self.wd,self.wp]

    def lift(self,y):
        return np.r_[np.repeat(y[:self.nold],self.repeats),y[self.nold:]]


def net(d,f,P,x,m,H,radial_order=8):
    chart,y=d.chart(x,m,radial_order)
    c,A,B=event(chart,f,d)
    energies=np.concatenate(d.energies(x,m))
    if abs(fsum(energies*c))>3e-13 or y>2e-14:
        raise ArithmeticError('Warm reference lost event-level energy')
    i=chart[3]
    return np.asarray([H*fsum(Q[i]*R) for Q in P for R in (A,B,A-B)]).reshape(len(P),3)


def run(c):
    x=c['x_end']
    base,y=unpack(c,False,x)
    if x not in (.01,.05):
        raise ArithmeticError('Use the declared early or late S4 benchmark endpoints')
    edges=base.ds.ed
    floor=float(min(edges[edges>=.01]))
    m,_,_,H=base.state(x,y)
    a=load('fine')[0]
    g=SPECS['fine']
    tab=[]
    def append(d,f,s,pressure):
        if not np.array_equal(d.wp,base.wp) or not np.array_equal(d.wd,base.wd) or (
            abs(fsum(d.ds.w*f[:d.ns])/fsum(base.ds.w*y[:base.ns])-1)>5e-14):
            raise ArithmeticError('Warm quadrature changed the charged chart or scalar number')
        P=[np.ones(d.ns),pressure]
        v=net(d,f,P,x,m,H,8)
        low=net(d,f,P,x,m,H,4)
        err=float(np.max(abs(v-low)/np.maximum(abs(v),1e-95)))
        if err>2e-6:
            raise ArithmeticError('The charge radial quadrature failed to resolve')
        tab.append(dict(rule=s,scalar_nodes=d.ns,number_and_warm_pressure_affine=v.tolist(),
                        charge_radial_mesh_drift=err))
    append(base,y,'geometric_W1',cell_weights(edges,m,floor))
    for j in (2,4):
        d=ScalarSubshellChart(a,g,2,4,1e-5,warm_floor=.01,warm_div=j)
        f=np.r_[y[np.searchsorted(edges,d.ds.p)-1],y[base.ns:]]
        append(d,f,f'geometric_W{j}',cell_weights(d.ds.ed,m,floor))
    d=ScalarGauss2(a,g,2,4,1e-5,gauss_floor=.01)
    q=d.ds.p
    append(d,d.transfer_step(y),'positive_Gauss2_volume',np.where(q>=floor,q*q/np.hypot(q,m),0.))
    for j in (8,16,32,64):
        d=WarmVolume(a,g,j)
        q=d.ds.p
        append(d,d.lift(y),f'Legendre_volume_{j}',np.where(q>=floor,q*q/np.hypot(q,m),0.))
    vals=[np.array(h['number_and_warm_pressure_affine']) for h in tab]
    T=vals[-1]
    err=[float(np.max(abs(a/b-1))) for a,b in zip(vals[-3:-1],vals[-2:],strict=True)]
    if max(err)>4e-5 or min(T[1,:2])<=0:
        raise ArithmeticError(('Higher warm volume integral still needs refinement',err))
    for d,j in zip(tab,vals,strict=True):
        d['relative_net_warm_pressure_vs_volume64']=float(j[1,2]/T[1,2]-1)
    return dict(x=x,complete_warm_cell_floor=floor,mass_held=m,clock_held=H,
                ref_16_to_32_and_32_to_64_relative_defects=err,results=tab,
                scope='Frozen constant occupations; warm pressure collision source only, no time evolution or physical pressure.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('S4',type=Path)
    print(json.dumps(run(read(p.parse_args().S4)),indent=2,allow_nan=False))
