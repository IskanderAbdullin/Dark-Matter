"""Energy- and species-number-conserving joint fermion momentum comparator."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import root
from scipy.special import expit

from early_pole_scalar_chemical_sensitivity_reference import Chemical


class Joint(Chemical):
    def __init__(self,shape,ratio,mix,kind):
        super().__init__(*shape,ratio,0.)
        self.mix=mix
        self.kind=kind

    def target(self,x,b,p,jac=False):
        nodes=np.r_[self.daughter_energies(x),self.parent_energies(x)]
        wb=np.r_[self.d.wf,self.pw]
        E=np.asarray(nodes)
        U=np.column_stack((np.arange(len(wb))<self.nf,
                           np.arange(len(wb))>=self.nf,E))
        B=(wb[:,None]*U).T
        state=np.r_[b,p]
        moments=B@state
        sc=np.maximum(abs(moments),1e-20)
        def f(z):
            return expit(U@z)
        def fd(z):
            q=U@z
            return U*(expit(q)*expit(-q))[:,None]
        nb,nn,_=moments
        guess=(np.log(nb/.365),np.log(nn/.365),-1.)
        sol=root(lambda z:(B@f(z)-moments)/sc,guess,
                 jac=lambda z:(B@fd(z))/sc[:,None],tol=8e-12)
        f0=f(sol.x)
        if np.max(abs((B@f0-moments)/sc))>5e-12:
            raise RuntimeError((sol.message,moments.tolist(),sol.x.tolist()))
        if jac:
            G=fd(sol.x)
            return G@np.linalg.solve(B@G,B)
        return f0,sol.x

    def at(self,x,m,y,jac=False):
        v=super().at(x,m,y,jac)
        if not self.mix or self.kind=='none':
            return v
        if self.kind=='parent':
            return self.apply_parent(x,y,v,jac)
        b=y[self.ns:self.iP];p=y[self.iP:]
        if jac:
            J=self.target(x,b,p,True)
            v[self.ns:,self.ns:]+=self.mix*(J-np.eye(self.nf+self.np))
        else:
            f,meta=self.target(x,b,p)
            v[self.ns:]+=self.mix*(f-y[self.ns:])
        return v

    def apply_parent(self,x,y,v,jac):
        p=y[self.iP:]
        if jac:
            q=self.elastic_target(x,p,True)-np.eye(self.np)
            v[self.iP:,self.iP:]+=self.mix*q
        else:
            v[self.iP:]+=self.mix*(self.elastic_target(x,p)-p)
        return v

    def starting_checks(self,s):
        x=s['x_end']
        b=np.asarray(s['endpoint_fermion_occupations'])
        p=np.asarray(s['endpoint_parent_occupations'])
        v=np.r_[b,p]
        f,fit=self.target(x,b,p)
        D=self.target(x,b,p,True)
        e=np.r_[self.daughter_energies(x),self.parent_energies(x)]
        w=np.r_[self.d.wf,self.pw]
        B=np.array([w*(np.arange(len(w))<self.nf),w*(np.arange(len(w))>=self.nf),e*w])
        inv=(B@(f-v))/(B@v)
        direc=v*np.cos(np.arange(len(v)))
        dz=3e-5
        out=[]
        for sign in (-1,1):
            z=v+sign*dz*direc
            out.append(self.target(x,z[:self.nf],z[self.nf:])[0])
        nd=(out[1]-out[0])/(2*dz)
        er=max(abs(nd-D@direc))/max(abs(nd))
        assert er<2e-7 and max(abs(inv))<7e-13
        return dict(initial_joint_target_temperature_over_T=float(-1/fit[-1]),
                    initial_species_and_heat_residuals=inv.tolist(),
                    initial_joint_target_derivative_error=float(er))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--rate',type=float,default=100.)
    p.add_argument('--kind',choices=['none','parent','joint'],default='joint')
    p.add_argument('--xend',type=float,default=.052)
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    s=json.loads(raw)
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}[tuple(s['cells'])]
    d=Joint(shape,s['mass_ratio'],a.rate,a.kind)
    meta=d.starting_checks(s)
    z=d.from_snapshot(s,a.xend)
    z.update(meta,fermion_momentum_model=a.kind,mixing_rate_over_radiation_H=a.rate)
    print(json.dumps(z,indent=2,allow_nan=False))
