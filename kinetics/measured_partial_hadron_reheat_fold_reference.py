"""Measured vacuum spectral sub-bands on the prescribed high-reheat clock."""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import kve, roots_legendre

from conditional_relic_normalization_reference import conversion, proxy
from kloe_pion_band_portal_reference import POINT_MUON, table as pion_data
from scalar_count_reheat_target_reference import Fold
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import M, m, mp
from timelike_portal_onsager_coverage_reference import V, alpha, factor, width


def bes_data(path):
    # pdftotext -layout of arXiv v2, Table II. R errors are stat and syst.
    pat=r"(?m)^\s*(2\.\d{4}|3\.\d{4})\s+(?:[\d.]+\s+){6}(2\.\d{3})\s*±\s*(0\.\d{3})\s*±\s*(0\.\d{3})\s*$"
    a=np.asarray(re.findall(pat,Path(path).read_text()),float).T
    assert a.shape==(4,14) and np.all(np.diff(a[0])>0.)
    assert np.array_equal(a[:2,[0,-1]],[[2.2324,3.6710],[2.286,2.405]])
    return a


def spectra(kloe,bes,n,mode):
    x,w=roots_legendre(n)
    s,cs,err=pion_data(kloe)
    q=s[:,None]+.005*x
    A=np.broadcast_to(.005*w,q.shape).ravel()
    k=(q.ravel(), A*np.repeat(cs,n)*q.ravel()/POINT_MUON,
                   A*np.repeat(err,n)*q.ravel()/POINT_MUON)
    a,r,e1,e2=bes_data(bes)
    groups=[]
    # Never interpolate through the J/psi strip between 3.08 and 3.4 GeV.
    for sl in (slice(None,6),slice(6,None)):
        E=a[sl]; c=r[sl]; e=(e1+e2)[sl]
        lo,hi=E[:-1,None]**2,E[1:,None]**2
        s=((hi+lo)+(hi-lo)*x)/2
        v=np.sqrt(s)
        if mode=="linear":
            Y=np.interp(v,E,c); Z=np.interp(v,E,e)
        else:
            edge=(E[:-1]+E[1:])/2
            # Exact split at the midpoint of adjacent energy points.
            low=np.c_[E[:-1],edge].ravel()[:,None]**2
            high=np.c_[edge,E[1:]].ravel()[:,None]**2
            s=((high+low)+(high-low)*x)/2
            lo,hi=low,high
            Y=np.repeat(np.c_[c[:-1],c[1:]].ravel()[:,None],n,axis=1)
            Z=np.repeat(np.c_[e[:-1],e[1:]].ravel()[:,None],n,axis=1)
        u=(hi-lo)/2*w
        groups.append((s.ravel(),(u*Y).ravel(),(u*Z).ravel()))
    h=tuple(np.concatenate([z[i] for z in groups]) for i in range(3))
    # Narrow J/psi sensitivity, rounded PDG mass, width and branching.
    Em=3.096900
    xi=9*np.pi*Em*(92.6e-6*.05971)*.877/(alpha*alpha)
    jj=(np.array([Em**2]),np.array([xi]),np.array([0.]))
    l,r=np.array([9.986,10.006])**2
    sl=(r+l)/2+(r-l)/2*x
    dw=(r-l)/2*w
    # One measured bracket anchor; a flat window is only a proxy.
    cleo=(sl,dw*3.497,dw*(.004+.064+.043))
    return (k,h,jj,cleo)


def integrands(T,packs):
    result=[]
    for s,R,err in packs:
        E=np.sqrt(s)
        K1=kve(1,E/T)*np.exp(-E/T)
        K0=kve(0,E/T)*np.exp(-E/T)
        z=[]
        for mass in (m,M):
            W=factor*(s+2*mass**2)*E**3*np.sqrt(1-4*mass**2/s)
            W/=((s-V**2)**2+(V*width)**2)
            g=W*T*K1
            en=W*(2*T*T*K1+T*E*K0)
            z.append([2*(R@g),R@en,2*(err@g),err@en])
        result.append(np.asarray(z))
    return np.asarray(result)


class HadronFold:
    def __init__(self,base,packs,ot):
        self.base=base; self.packs=packs
        self.x,self.w=roots_legendre(ot)
        self.edges=self.base.edges
        y=[np.zeros((len(self.packs),2,4))]
        for a,b in zip(self.edges[:-1],self.edges[1:]):
            y.append(y[-1]+self.panel(a,b))
        self.cum=np.asarray(y)

    def panel(self,a,b):
        z=np.zeros((len(self.packs),2,4))
        for u,w in zip((a+b)/2+(b-a)/2*self.x,wid(self.w,a,b)):
            T=np.exp(u)
            s,rho=self.base.densities(u)
            H=np.sqrt(rho/3)/MPL
            dt=w*(3+self.base.bg.sg.derivative()(u))/(3*H)
            G=integrands(T,self.packs)
            z+=G*np.array([dt/s,dt*(self.base.sf/s)**(4/3)/self.base.rf,
                           dt/s,dt*(self.base.sf/s)**(4/3)/self.base.rf])
        return z

    def sample(self,T):
        u=np.log(T); j=np.searchsorted(self.edges,u,side="right")-1
        return self.base.sample(T),self.cum[j]+self.panel(self.edges[j],u)


def wid(w,a,b):
    return (b-a)/2*w


def report(eos,kloe,bes):
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    base=Fold(eos,(8,21,16),6)
    controls=[]
    for n,mode in ((8,"linear"),(14,"piecewise")):
        F=HadronFold(base,spectra(kloe,bes,n,mode),6)
        memo={}
        def val(T,case,error=0.):
            if T not in memo:memo[T]=F.sample(T)
            leps,meas=memo[T]
            L=leps[:2]+(meas[:case,:,0]+error*meas[:case,:,2]).sum(axis=0)
            W=leps[2:].sum()+(meas[:case,:,1]+error*meas[:case,:,3]).sum()
            o=proxy([L[1],L.sum()],f)
            return o[0]+o[1]*mp/m,W,o
        rows=[]
        for c,name in enumerate(("leptons","plus KLOE pion band","plus BESIII measured continuum sub-bands","plus illustrative narrow J/psi","plus conditional near-pole CLEO point proxy")):
            for frac in (.5,.99,1.):
                T=brentq(lambda t: val(t,c)[0]-frac*.12,.2,.6,xtol=1e-8)
                def br(err):return brentq(lambda t:val(t,c,err)[0]-frac*.12,.2,.6,xtol=1e-8)
                _,e,v=val(T,c)
                leps,h=memo[T]
                channels=mp*f*(h[:,0,0]+2*h[:,1,0])
                rows.append(dict(convention=name,target_mass_fraction=frac,T_R_GeV=T,
                    coherent_input_error_shift_bounds_GeV=[br(1.),br(-1.)],
                    stable_no_loss_Omega_h2=v[1],scalar_direct_prompt_Omega_h2=v[0],
                    partial_hadron_count_contributions_omega_if_all_convert=channels[:c].tolist(),
                    energy_floor_over_radiation_at_1MeV=e))
        controls.append(dict(interpolation=mode,invariant_nodes_per_interval=n,rows=rows))
    # KLOE published local MB check at 0.3 GeV, using the same invariant normalization.
    from kloe_pion_band_portal_reference import row
    ours=integrands(.3,spectra(kloe,bes,8,"linear"))[0,0]
    K=row(.3,*pion_data(kloe))["measured_band_N_Q_Q2"]
    num=max(abs(ours[0]/(2*K[0])-1),abs(ours[1]/K[1]-1))
    assert num<2e-10
    from combined_timelike_spacelike_fd_modes_reference import eval_modes,chi_d,PositivePortal,evaluate,StandardModel
    p=PositivePortal((10,28,20)); E=StandardModel(eos)
    loc=[]
    for q in controls[0]["rows"]:
        if q["target_mass_fraction"]!=.99: continue
        T=q["T_R_GeV"]
        d=evaluate(T,0.)
        v=eval_modes(T,d["rho"],chi_d(T,d),E,p.moments(T,T,0.)[3:6])
        loc.append(dict(convention=q["convention"],T_GeV=T,
                        lepton_only_local_FD_mode_rates_over_H=v.tolist()))
    return dict(controls=controls,independent_KLOE_kernel_relative_error=float(num),
                lepton_only_local_thermal_checks_at_99_percent_thresholds=loc)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("eos","kloe_combined_data","besiii_arxiv_pdftotext_layout"):p.add_argument(k)
    a=p.parse_args()
    out=report(a.eos,a.kloe_combined_data,a.besiii_arxiv_pdftotext_layout)
    out["data_hashes"]={Path(k).name:hashlib.sha256(Path(k).read_bytes()).hexdigest()
                       for k in (a.kloe_combined_data,a.besiii_arxiv_pdftotext_layout)}
    out["scope"]="Only the specified disjoint vacuum hadronic sub-bands plus a separate isolated J/psi convention. An optional CLEO point proxy is not a measured continuous band. They are added with classical thermal weights to the FD lepton source on a prescribed radiation background. No complete hadronic spectral function, finite-temperature QCD, inverse depletion or self-consistent relic inference."
    print(json.dumps(out,indent=2,allow_nan=False))
