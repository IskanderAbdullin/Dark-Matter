"""Exploratory frozen cold-branch flux with a large-N-inspired vertex.

Neither the N=1 real-scalar collision kernel nor a full plasma damping rate.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_scalar_elastic_stencil_reference import Elastic
from early_pole_scalar_nr_vertex_response_reference import primitive


class Screen:
    def __init__(self,edges,f,m,lam,order):
        self.edges=np.asarray(edges)
        self.f=np.asarray(f)
        self.m=m
        self.g=lam/(8*m*m)
        self.angle,self.w=roots_legendre(order)
        self.azimuth=np.cos(np.pi*(self.angle+1))

    def vertex(self,P,omega):
        P=np.maximum(P,2e-13)
        z=2*self.m*omega/P**2
        l=self.edges[:-1][None,:]/P[:,None]
        h=self.edges[1:][None,:]/P[:,None]
        a=np.c_[(-1-z)/2,(-1+z)/2,(1-z)/2,(1+z)/2]
        c=np.array([1.,1.,-1.,-1.])
        part=primitive(h[:,:,None],a[:,None,:])-primitive(l[:,:,None],a[:,None,:])
        re=np.sum((part@c)*self.f,axis=1)
        r1=np.minimum(abs(1-z),abs(1+z))[:,None]*.5
        r2=np.maximum(abs(1-z),abs(1+z))[:,None]*.5
        lo=np.maximum(l,r1)
        hi=np.minimum(h,r2)
        im=np.pi*np.sum(.5*np.maximum(hi-lo,0.)*(hi+lo)*self.f,axis=1)*np.sign(z)
        norm=self.m*self.g*P/(4*np.pi**2)
        return 1./((1+norm*re)**2+(norm*im)**2)

    def __call__(self,p,k,u,v,delta,s,lab,e1,e2,outK):
        a=np.maximum(lab,1e-30)
        z1=(p*p+p*k*u)/a
        t1=p*k*np.sqrt(1-u*u)/a
        pc=np.sqrt(delta/2)
        z3=.5*lab-(e1+e2)/np.sqrt(s)*pc*v
        t3=pc*np.sqrt(1-v*v)
        outq2=outK*(2*self.m+outK)
        omega=(p*p-outq2)/(2*self.m)
        normal=(z1-z3)**2+(t1-t3)**2
        mixed=2*t1*t3
        acc=np.zeros(len(p))
        for cos,w in zip(self.azimuth,self.w,strict=True):
            transfer=np.sqrt(np.maximum(normal+mixed*(1-cos),0.))
            acc+=.5*w*self.vertex(transfer,omega)
        return acc


def one(s,order,azimuth):
    m=s['checkpoints'][-1]['scalar_pole_over_T']
    edges=np.array(s['scalar_momentum_edges'])
    q=np.array(s['scalar_momentum_nodes'])
    F=np.maximum(np.array(s['endpoint_scalar_occupations']),0.)
    r=np.searchsorted(edges,1e-5)
    assert edges[r]==1e-5
    q=q[:r]; edges=edges[:r+1]; F=F[:r]
    lm=7.65175e-8
    screen=Screen(edges,F,m,lm,azimuth)
    # Compare exactly the same closed IR grid and reversal graph.
    A=Elastic(q,m,edges,order,lm)
    B=Elastic(q,m,edges,order,lm,screen)
    assert A.weight.shape==B.weight.shape
    u,ha=A.flux(F);v,hb=B.flux(F)
    c=A.vol@F
    h=np.sqrt(np.pi**2*106.75/90)*(86./s['x_end'])/2.435e18
    u0=np.sum(abs(u));v0=np.sum(abs(v))
    assert abs(v.sum())<3e-12*v0
    assert abs(v@B.K)<3e-12*(abs(v)@B.K)
    assert ha>0 and hb>0
    return dict(x=s['x_end'],mass_ratio=s['mass_ratio'],IR_cells=len(q),angular_order=order,azimuth_order=azimuth,
                IR_projected_events=A.weight.size,raw_geometry_out_of_grid=A.raw_escape,
                contact_g_over_Tref_minus2=screen.g,
                mean_vertex_sq_ratio_over_represented_dilute_weights=float(B.weight.sum()/A.weight.sum()),
                formal_bare_IR_redistribution_rate_over_Href=float(u0/(c*h)),
                formal_screened_cold_redistribution_rate_over_Href=float(v0/(c*h)),
                formal_net_IR_shape_flux_ratio=float(v0/u0),
                formal_entropy_production_ratio=float(hb/ha))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    args=p.parse_args()
    out={}
    for name in args.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d=json.loads(raw)
        out[f.name]=[one(d,angle,az) for angle,az in [(6,6),(10,10)]]
    print(json.dumps(out,indent=2,allow_nan=False))
