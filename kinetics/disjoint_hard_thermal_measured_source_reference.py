"""Prescribed hybrid source with no double count of the measured mass windows."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

from born_thermal_quark_hard_portal_reference import (
    XEDGES,UEDGES,gauss,kernel,local,thermal_trace_over_vacuum)
from conditional_relic_normalization_reference import conversion
from measured_partial_hadron_reheat_fold_reference import Fold,HadronFold,spectra
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp

# KLOE support is entirely below x=4 for these T. BESIII uses two
# separately measured intervals, without interpolation through J/psi.
WINDOWS=np.array([[2.2324,3.08],[3.4,3.6710]])


class HardComplement:
    def __init__(self,base,tmin=.25,tmax=.48,n=6,x_edges=None):
        self.base,self.n=base,n
        self.a,self.b=np.log([tmin,tmax])
        self.mass_edges=np.array(XEDGES if x_edges is None else x_edges,dtype=float)
        if (self.mass_edges.ndim!=1 or len(self.mass_edges)<2 or self.mass_edges[0]!=4 or
                not np.all(np.isfinite(self.mass_edges)) or not np.all(np.diff(self.mass_edges)>0)):
            raise ValueError("Born mass panels must increase from x=4.")
        cross=(WINDOWS[:,:,None]/self.mass_edges).ravel()
        self.edges=np.unique(np.clip(np.r_[self.a,self.b,base.edges,np.log(cross)],self.a,self.b))
        u,du=gauss(UEDGES,n)
        self.u,self.du=u,du
        self.cum=np.cumsum(np.r_[np.zeros((1,2,2)),[self.panel(a,b) for a,b in zip(self.edges[:-1],self.edges[1:])]],axis=0)

    def rates(self,T):
        q=WINDOWS/T
        edges=np.unique(np.clip(np.r_[self.mass_edges,q.ravel()],self.mass_edges[0],self.mass_edges[-1]))
        x,dx=gauss(edges,self.n)
        not_measured=~np.any((x[:,None]>q[:,0])&(x[:,None]<q[:,1]),axis=1)
        x,dx=x[:,None],dx[:,None]
        u=self.u
        e=x+u
        R=2*thermal_trace_over_vacuum(x,u)
        W=2*dx*self.du*np.sqrt(u*(2*x+u))/np.expm1(e)
        N=(R*W).sum(axis=1)
        Q=(R*W*e).sum(axis=1)
        kin=(x[:,0],np.stack([N,N*not_measured]),np.stack([Q,Q*not_measured]))
        return local(T,kin)

    def panel(self,a,b):
        if b<=a:
            return np.zeros((2,2))
        u,w=gauss([a,b],self.n)
        answer=np.zeros((2,2))
        for v,dv in zip(u,w):
            T=np.exp(v)
            s,r=self.base.densities(v)
            H=np.sqrt(r/3)/MPL
            dt=dv*(3+self.base.bg.sg.derivative()(v))/(3*H)
            shell=self.rates(T)
            answer[:,0]+=dt/s*(shell[:,0,0]+2*shell[:,1,0])
            answer[:,1]+=dt*(self.base.sf/s)**(4/3)/self.base.rf*shell[:,:,1].sum(axis=1)
        assert np.all(answer[0]>=answer[1])
        return answer

    def sample(self,tr):
        u=np.log(tr)
        if u<=self.a:
            return np.zeros((2,2))
        assert u<self.b
        j=np.searchsorted(self.edges,u,side="right")-1
        return self.cum[j]+self.panel(self.edges[j],u)


def run(eos,kloe,bes):
    A=Fold(eos,(8,21,16),6,upper=.5)
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]*mp
    V=HadronFold(A,spectra(kloe,bes,10,"linear"),6)
    def measured(T):
        lept,h=V.sample(T)
        q=lept[:2]+h[:2,:,0].sum(axis=0)
        return np.array([f*(q[0]+2*q[1]),lept[2:].sum()+h[:2,:,1].sum()])
    root=brentq(lambda T:measured(T)[0]-.1188,.28,.4,xtol=2e-9)
    results=[]
    independent=[]
    for tmin in (.25,.30):
        refinements=[]
        for n in (4,8):
            H=HardComplement(A,tmin,n=n)
            if n==8:
                K=kernel(8)
                for T in (.25,.30,.35,.42):
                    R=local(T,K)[0]
                    independent.append(float(np.max(np.abs(H.rates(T)[0]/R-1.))))
            def h(T):
                a=H.sample(T).copy()
                a[:,0]*=f
                return a
            root2=brentq(lambda T:measured(T)[0]+h(T)[1,0]-.1188,.27,.40,xtol=2e-9)
            full,comp=h(root)
            val=measured(root2)+h(root2)[1]
            refinements.append(dict(nodes_per_panel=n,hybrid_T_R_for_99_percent_optimistic_count_GeV=root2,
                 at_measured_only_root=dict(full_Born=full.tolist(),nonoverlapping_Born=comp.tolist(),
                                            fractional_Born_overlap=((full-comp)/full).tolist()),
                 hybrid_energy_floor_at_its_root=val[1]))
        assert abs(refinements[0]["hybrid_T_R_for_99_percent_optimistic_count_GeV"]/
                   refinements[1]["hybrid_T_R_for_99_percent_optimistic_count_GeV"]-1)<2e-5
        results.append(dict(Born_lower_visible_T_GeV=tmin,refinements=refinements))
    assert 4*.25>.975 and max(independent)<3e-9
    return dict(measured_only_99_percent_optimistic_root_GeV=root,
         excluded_BESIII_M_GeV=WINDOWS.tolist(),rows=results,
         max_full_Born_kernel_repartition_relative_defect=max(independent),
         input_sha256={Path(p).name:hashlib.sha256(Path(p).read_bytes()).hexdigest()
                       for p in (eos,kloe,bes)},
         columns_in_hard_arrays=["optimistic all-converted scalar Omega h^2","dark energy floor/visible at 1 MeV"],
         scope="Conditional patchwork: lepton FD source plus specified KLOE/BESIII vacuum bands, and massless thermal Born only away from the measured invariant-mass windows over bounded x,u and T. No overlap counting, interacting thermal QCD, near-pole proxy, charmonium resonance input, inverse or evolution. Optimistic count, not a physical relic threshold or bound on exact QCD.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for key in ("eos","kloe","besiii_pdftotext_layout"):
        p.add_argument(key)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.besiii_pdftotext_layout),indent=2,allow_nan=False))
