"""Conditional strict-NLO increment on measured perturbative strips, not full thermal QCD."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator,interp1d
from scipy.optimize import brentq

from born_thermal_quark_hard_portal_reference import XEDGES,UEDGES,gauss,local
from conditional_relic_normalization_reference import conversion
from disjoint_hard_thermal_measured_source_reference import HardComplement,WINDOWS
from hard_vector_nlo_strip_reference import spectra
from measured_partial_hadron_reheat_fold_reference import Fold,HadronFold,spectra as measured_spectra
from three_species_decay_onsager_reference import mp


def table(path):
    a=np.atleast_2d(np.loadtxt(path))
    X=np.unique(a[:,0])
    if X.size<3 or X[0]!=4. or X[-1]>12. or np.max(np.diff(X))>1.:
        raise ValueError("Use complete contiguous mass strips from x=4, with separation at most one and x<=12.")
    if not np.all(np.isfinite(a)) or not np.all(np.diff(a[:,0])>=0.):
        raise ValueError("Nonfinite or unordered native input.")
    assert len(a)==len(X)*len(UEDGES)
    pp=a[:,1].reshape(len(X),-1)
    V=3*X[:,None]**2/(4*np.pi)
    u,lo=spectra(X[:,None],pp)
    assert np.max(abs(u-UEDGES))<1e-12
    assert np.max(abs(lo/a[:,2].reshape(u.shape)-1.))<1e-11
    return X,np.asarray(UEDGES),a[:,3].reshape(u.shape)/V


def interp(X,U,C,mode):
    cls=PchipInterpolator if mode=="pchip" else interp1d
    alongU=cls(U,C,axis=1)
    return lambda u:cls(X,alongU(u),axis=0)


class Correction(HardComplement):
    def __init__(self,base,table_values,tmin,n,mode):
        x,u,c=table_values
        un,_=gauss(UEDGES,n)
        self.cf=interp(x,u,c,mode)(un)
        self.mass_edges=x
        self.base,self.n=base,n
        self.a,self.b=np.log([tmin,.48])
        cross=(WINDOWS[:,:,None]/np.union1d(XEDGES,x)).ravel()
        self.edges=np.unique(np.clip(np.r_[self.a,self.b,base.edges,np.log(cross)],self.a,self.b))
        self.u,self.du=gauss(UEDGES,n)
        self.cum=np.cumsum(np.r_[np.zeros((1,2,2)),
            [self.panel(a,b) for a,b in zip(self.edges[:-1],self.edges[1:])]],axis=0)

    def rates(self,T):
        q=WINDOWS/T
        edges=np.unique(np.clip(np.r_[self.mass_edges,q.ravel()],self.mass_edges[0],self.mass_edges[-1]))
        x,dx=gauss(edges,self.n)
        outside=~np.any((x[:,None]>q[:,0])&(x[:,None]<q[:,1]),axis=1)
        u=self.u
        C=self.cf(x)
        e=x[:,None]+u
        # For three unit-color quark flavors with charge factor sum=2,
        # C is the strict NLO trace relative to the one-flavor vacuum
        # trace, with g_s^2 factored out until the final report.
        R=2*C
        W=2*dx[:,None]*self.du*np.sqrt(u*(2*x[:,None]+u))/np.expm1(e)
        N=(R*W).sum(axis=1)
        Q=np.sum(R*W*e,axis=1)
        assert np.all(N>0.)
        return local(T,(x,np.stack([N,N*outside]),np.stack([Q,Q*outside])))


def mid_checks(table_values,path):
    if path is None:
        return []
    X,U,C=table_values
    rows=[]
    for x,k,lo,c in np.atleast_2d(np.loadtxt(path)):
        u,l=spectra(x,k)
        assert abs(lo/l-1)<5e-11 and X[0]<x<X[-1]
        t=c/(3*x*x/(4*np.pi))
        d={}
        for mode in ("linear","pchip"):
            y=interp(X,U,C,mode)(np.array([u]))(np.array([x]))[0,0]
            d[mode]=float(y/t-1.)
        rows.append(dict(x=x,u=float(u),prediction_relative_errors=d))
    return rows


def run(eos,kloe,bes,probe,check):
    K=table(probe)
    F=Fold(eos,(8,21,16),6,upper=.50)
    V=HadronFold(F,measured_spectra(kloe,bes,10,"linear"),6)
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    target=.1188
    def measured(T):
        l,h=V.sample(T)
        q=l[:2]+h[:2,:,0].sum(axis=0)
        return np.array([f*(q[0]+2*q[1]),l[2:].sum()+h[:2,:,1].sum()])
    results=[]
    for tmin in (.25,.30):
        B=HardComplement(F,tmin=tmin,n=8)
        lo=brentq(lambda T:measured(T)[0]+f*B.sample(T)[1,0]-target,.27,.40,xtol=2e-9)
        rows=[]
        for mode,n in (("linear",4),("pchip",8)):
            C=Correction(F,K,tmin,n,mode)
            vals=[]
            for a in (.2,.3,.4):
                g2=4*np.pi*a
                def rate(T):
                    z=B.sample(T)[1]+g2*C.sample(T)[1]
                    return measured(T)+z*np.array([f,1.])
                hi=brentq(lambda T:rate(T)[0]-target,.27,.40,xtol=2e-9)
                c=C.sample(lo)*g2*np.array([f,1.])
                vals.append(dict(alpha_s=a,T_R_GeV=hi,root_shift_GeV=hi-lo,
                    incremental_omega_at_baseline_without_overlap=float(c[1,0]),
                    incremental_energy_floor_at_baseline_without_overlap=float(c[1,1]),
                    fraction_of_strict_increment_cut_for_overlap=float(1-c[1,0]/c[0,0])))
            rows.append(dict(interpolation=mode,nodes_per_panel=n,alpha_cases=vals))
        results.append(dict(Born_T_min_GeV=tmin,measured_plus_disjoint_Born_root_GeV=lo,results=rows))
    return dict(controls=results,mass_over_T_interval=K[0][[0,-1]].tolist(),
                independent_native_point_interpolation_checks=mid_checks(K,check),
                data_hashes={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,kloe,bes,probe,check) if p is not None},
                scope="Conditional addition to the previous disjoint hard Born comparator, using strict NLO only in the supplied bounded mass interval with resolved boost probes and finite-band interpolation. Same radiation clock, same optimistic conversion. Not a complete finite-T QCD spectrum, QCD error bars, a full NLO root, a physical dark-matter abundance, nor a production fit.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","kloe","besiii","complete_strips"):
        p.add_argument(k)
    p.add_argument("--intermediate-check")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.besiii,a.complete_strips,a.intermediate_check),indent=2,allow_nan=False))
