"""Reversible cold 4-to-2 threshold stencil at fixed mass, one reaction family."""
import argparse
import gzip
from collections import Counter
from itertools import combinations_with_replacement
import json
from math import factorial, prod
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

from early_pole_scalar_number_change_threshold_reference import LAMBDA, threshold


class NumberChange:
    def __init__(self,q,edges,m,qlo,amp_scale):
        self.E=np.hypot(q,m)
        vol=np.diff(edges**3)/(6*np.pi**2)
        cold=np.flatnonzero((edges[:-1]>=qlo)&(edges[1:]<=1e-5))
        I=np.array(list(combinations_with_replacement(cold,4)))
        targets=self.E[I].sum(axis=1)/2
        k=np.searchsorted(self.E,targets,side='right')
        assert np.min(k)>0 and np.max(k)<len(self.E)
        r=(targets-self.E[k-1])/(self.E[k]-self.E[k-1])
        full=np.c_[I,k-1,k]
        signed=np.c_[-np.ones_like(I),2*(1-r),2*r]
        active=np.unique(full)
        self.active=active
        jj=np.searchsorted(active,full)
        rows=np.repeat(np.arange(len(I)),6)
        shape=(len(I),len(active))
        self.S=coo_matrix((signed.ravel(),(rows,jj.ravel())),shape=shape).tocsr()
        self.A=coo_matrix((abs(signed).ravel(),(rows,jj.ravel())),shape=shape).tocsr()
        self.E=self.E[active]
        self.vol=vol[active]
        self.cold=np.isin(active,cold)
        self.W=amp_scale*LAMBDA**4/m**8*threshold()['event_coefficient_times_lambda4_over_m8']
        multiplicity=np.array([factorial(4)/prod(factorial(x) for x in Counter(row).values()) for row in I])
        self.W*=np.prod(vol[I],axis=1)*multiplicity
        self.incoming=jj[:,:4]
        self.outgoing=jj[:,4:]
        self.outw=signed[:,4:]
        assert np.max(abs(self.S@self.E))<2e-14*m
        assert np.max(abs(self.S@np.ones(len(active))+2))<5e-15

    def events(self,f):
        f=np.maximum(f,0.)
        log=np.log(np.where(f>0,f,1.))
        h=np.logaddexp(0.,-log)
        D=self.S@h
        left=h[self.incoming].sum(axis=1)
        logamp=self.A@log+left+np.maximum(D,0.)
        R=np.exp(logamp)*np.sign(D)*(-np.expm1(-abs(D)))
        if np.any(f==0.):
            vin=(f[self.incoming]>0).all(axis=1)
            vout=((f[self.outgoing]>0)|(self.outw==0.)).all(axis=1)
            forward=vin&(~vout); reverse=(~vin)&vout
            fc=f[self.outgoing]
            left_f=np.log1p(f[self.incoming]).sum(axis=1)
            right_f=(self.outw*np.log1p(fc)).sum(axis=1)
            h0=self.A@log
            R[forward]=np.exp((h0+right_f)[forward])
            R[reverse]=-np.exp((h0+left_f)[reverse])
            R[~vin&~vout]=0
            D[forward]=np.inf
            D[reverse]=-np.inf
            D[~vin&~vout]=0.
        return R*self.W,D

    def flux(self,f):
        R,D=self.events(f)
        positive=(R!=0)
        assert R[positive]@D[positive]>=0
        return self.S.T@R

    def verify(self,m):
        G=self.vol
        th=.7*m
        thermal=1./np.expm1(self.E/th)
        above=1./np.expm1((self.E-.13*m)/th)
        below=1./np.expm1((self.E+.13*m)/th)
        ve=self.flux(thermal)
        vc=self.flux(above)
        vo=self.flux(below)
        R,_=self.events(above)
        incl=abs(vc.sum()+2*R.sum())/max(np.sum(abs(vc)),1e-300)
        therm=np.sum(abs(ve))/np.sum(abs(vc))
        assert therm<1e-12 and incl<1e-12
        assert vc.sum()<0 and vo.sum()>0
        trial=above.copy()
        ix=np.flatnonzero(~self.cold)[0]
        trial[ix]=0.
        assert self.flux(trial)[ix]>0.
        return dict(thermal_residual_relative=float(therm),stoichiometry_residual=float(incl),
                    positive_chemical_potential_removes_particles=True)


def run(spec,qlo):
    m=spec['checkpoints'][-1]['scalar_pole_over_T']
    edges=np.array(spec['scalar_momentum_edges'])
    q=np.array(spec['scalar_momentum_nodes'])
    f=np.maximum(np.array(spec['endpoint_scalar_occupations']),0.)
    C=NumberChange(q,edges,m,qlo,1.)
    checks=C.verify(m)
    n=C.vol*f[C.active]
    N=n.sum()
    f0=n/C.vol
    v0=C.flux(f0)
    norm=-v0.sum()/N
    assert norm>0 and np.min(n)>0
    y0=n/N
    def rhs(_,y):
        return C.flux(np.maximum(y,0.)*N/C.vol)/(norm*N)
    y=solve_ivp(rhs,(0.,.03),y0,t_eval=np.r_[0.,np.geomspace(1e-10,.03,45)],
                method='BDF',rtol=3e-9,atol=np.maximum(y0*2e-11,1e-18))
    assert y.success,y.message
    Y=y.y
    assert np.min(Y)>0
    residual=np.max(abs((C.E@Y)/(C.E@y0)-1))
    assert residual<2e-9
    F=Y*N/C.vol[:,None]
    S=C.vol@(np.log1p(F)+F*np.log1p(1/F))
    assert min(np.diff(S))>=-1e-12*S[-1]
    nall=C.vol@f0
    H=np.sqrt(np.pi**2*106.75/90)*(86./spec['x_end'])/2.435e18
    return dict(x=spec['x_end'],mass_ratio=spec['mass_ratio'],retained_incoming_lower_cut_over_T=qlo,
                cold_cells=int(C.cold.sum()),active_cells=len(C.cold),reversible_reactions=C.S.shape[0],
                initial_number_fraction_represented=float(nall/np.sum(f*np.diff(edges**3)/(6*np.pi**2))),
                initial_formal_number_loss_rate_over_Href=float(norm/H),
                dimensionless_time=[float(t) for t in y.t],
                active_number_fraction_remaining=Y.sum(axis=0).tolist(),
                active_cold_fraction_of_initial_total=Y[C.cold].sum(axis=0).tolist(),
                max_energy_drift=float(residual),entropy_relative_increase=float(S[-1]/S[0]-1),
                static=checks)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectrum')
    p.add_argument('--lower',type=float,default=1e-6)
    a=p.parse_args()
    path=Path(a.spectrum)
    raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
    print(json.dumps(run(json.loads(raw),a.lower),indent=2,allow_nan=False))
