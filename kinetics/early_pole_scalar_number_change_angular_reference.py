"""Reversible angular deposition for a cold threshold 4-to-2 amplitude."""
import argparse
import gzip
from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import qmc

from early_pole_scalar_number_change_stencil_reference import NumberChange
from early_pole_scalar_number_change_threshold_reference import LAMBDA, threshold


@lru_cache(maxsize=24)
def data(cold_key,q_key,power):
    cold=np.array(cold_key)
    q=np.array(q_key)
    I=np.array(list(combinations_with_replacement(cold,4)))
    mult=np.array([factorial(4)/prod(factorial(x) for x in Counter(row).values()) for row in I])
    u=qmc.Sobol(d=6,scramble=True,seed=50732).random_base2(power)
    c=2*u[:,[0,1,3]]-1.
    r=np.sqrt(1-c*c)
    a=np.column_stack((np.zeros(len(u)),2*np.pi*u[:,2],2*np.pi*u[:,4]))
    vec=np.stack((r*np.cos(a),r*np.sin(a),c),axis=2)
    dirs=np.concatenate((np.broadcast_to(np.array([0.,0.,1.]),(len(u),1,3)),vec),axis=1)
    lab=np.einsum('na,maj->nmj',q[I],dirs)
    dot=np.linalg.norm(lab,axis=2)
    return I,mult,dot,2*u[:,5]-1.


class AngularNumberChange(NumberChange):
    def __init__(self,q,edges,m,qlo,amp_scale,power=3):
        E=np.hypot(q,m)
        vol=np.diff(edges**3)/(6*np.pi**2)
        cold=np.flatnonzero((edges[:-1]>=qlo)&(edges[1:]<=1e-5))
        I,mult,lab,cos=data(tuple(cold),tuple(q),power)
        inc=E[I].sum(axis=1)
        shat=inc[:,None]**2-lab**2
        assert np.min(shat)>=(16-1e-10)*m*m
        phase=np.sqrt(np.maximum(1-4*m*m/shat,0.))
        shift=.5*lab*phase*cos
        high=(inc[:,None]/2+shift).ravel()
        low=(inc[:,None]/2-shift).ravel()
        triples=len(high)
        i=np.repeat(I,len(cos),axis=0)
        ix=np.searchsorted(E,low,side='right')
        iy=np.searchsorted(E,high,side='right')
        assert min(np.min(ix),np.min(iy))>0 and max(np.max(ix),np.max(iy))<len(E)
        wx=(low-E[ix-1])/(E[ix]-E[ix-1])
        wy=(high-E[iy-1])/(E[iy]-E[iy-1])
        support=np.c_[i,ix-1,ix,iy-1,iy]
        signed=np.c_[-np.ones_like(i),1-wx,wx,1-wy,wy]
        self.active=np.unique(support)
        jj=np.searchsorted(self.active,support)
        shp=(triples,len(self.active))
        rr=np.repeat(np.arange(triples),8)
        self.S=coo_matrix((signed.ravel(),(rr,jj.ravel())),shape=shp).tocsr()
        self.A=coo_matrix((abs(signed).ravel(),(rr,jj.ravel())),shape=shp).tocsr()
        self.E=E[self.active]
        self.vol=vol[self.active]
        self.cold=np.isin(self.active,cold)
        # The residual angle dependence of the tree amplitude is left at
        # threshold. Include exact on-shell phase space and initial energies.
        factor=(mult*np.prod(vol[I]*m/E[I],axis=1))[:,None]*phase/(len(cos)*np.sqrt(.75))
        self.W=factor.ravel()*amp_scale*LAMBDA**4/m**8*threshold()['event_coefficient_times_lambda4_over_m8']
        self.incoming=jj[:,:4]
        self.outgoing=jj[:,4:]
        self.outw=signed[:,4:]
        assert np.max(abs(self.S@self.E))<2e-14*m
        assert np.max(abs(self.S@np.ones(len(self.active))+2))<5e-15


def one(spec):
    r=spec['checkpoints'][-1]
    m=r['scalar_pole_over_T']
    q=np.array(spec['scalar_momentum_nodes'])
    edges=np.array(spec['scalar_momentum_edges'])
    f=np.maximum(np.array(spec['endpoint_scalar_occupations']),0.)
    rows=[]
    for order in (None,2,3,4,5):
        a=(NumberChange(q,edges,m,1e-6,1.) if order is None else
           AngularNumberChange(q,edges,m,1e-6,1.,order))
        V=a.flux(f[a.active])
        h=np.sqrt(np.pi**2*106.75/90)*(86/r['x'])/2.435e18
        nt=np.diff(edges**3)@f/(6*np.pi**2)
        raw=abs(V@a.E)/np.sum(abs(V)*a.E)
        assert raw<1e-12 and V.sum()<0
        rows.append(dict(angles=0 if order is None else 2**order,events=a.S.shape[0],
                         formal_number_loss_over_Href=float(-V.sum()/(nt*h)),
                         hot_receiving_nodes=q[a.active[~a.cold]].tolist(),
                         relative_energy_residual=float(raw),static=a.verify(m)))
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    a=p.parse_args()
    d={}
    for name in a.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d[f.name]=one(json.loads(raw))
    print(json.dumps(d,indent=2,allow_nan=False))
