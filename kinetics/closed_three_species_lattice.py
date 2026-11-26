"""CLOSED flat-space 3-species resonant-lattice test for chi1 <-> chi2 + phi.

Unlike fixed-parent-bath models, all three species evolve. A UNIFORM
daughter energy grid permits EXACT discrete reaction resonances:
 e_phi_i=(i+1/2)*Delta, e_chi_j=(j+1/2)*Delta,
 E_parent_{i+j}=(i+j+1)*Delta=e_phi_i+e_chi_j.

This makes BOTH independent particle-number invariants, the complete
discrete energy and a quantum H-theorem exact on the semidiscrete ODE.
The uniform grid underresolves IR Bose stimulation, and finite energy
cutoffs/representative centres approximate the continuum collision rate.
This is structural verification, NOT a physical cold-burst spectrum.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit


def allowed_area_from_zero(edges, M1_T):
    """Analytic area of a daughter-energy cell with q_phi*q_chi>=M1^2/4.

    Includes daughter energy cells starting at zero; no 1/0 poles occur
    because the partial-domain lower x is strictly positive.
    """
    edges=np.asarray(edges,dtype=float)
    if edges[0]<0 or not np.all(np.diff(edges)>0):
        raise ValueError("Need nonnegative strictly increasing energy edges")
    c=M1_T*M1_T/4
    xl=edges[:-1,None];xh=edges[1:,None]
    yl=edges[None,:-1];yh=edges[None,1:]
    limit=np.divide(c,yl,out=np.full_like(yl,np.inf),where=yl>0)
    xa=np.maximum(xl,c/yh)
    xb=np.minimum(xh,limit)
    part=np.where(xb>xa,yh*(xb-xa)-c*np.log(xb/xa),0.)
    full=np.maximum(xh-np.maximum(xl,limit),0.)*(yh-yl)
    return np.maximum(part+full,0.)


class ClosedResonantDecay:
    def __init__(self,nq=50,delta=0.08,M1_T=.1,g1=4,gchi=4):
        if nq<5 or delta<=0 or M1_T<=0:
            raise ValueError("Invalid resonant lattice")
        self.nq=nq;self.delta=delta;self.M1_T=M1_T
        self.edges=delta*np.arange(nq+1)
        self.e=(np.arange(nq)+.5)*delta
        self.vphi=(self.edges[1:]**3-self.edges[:-1]**3)/(6*math.pi**2)
        self.vchi=gchi*self.vphi
        self.parent_E=(np.arange(2*nq-1)+1)*delta
        elo=np.maximum((np.arange(2*nq-1)+.5)*delta,M1_T)
        ehi=np.maximum((np.arange(2*nq-1)+1.5)*delta,M1_T)
        # Mother Fermi sphere phase-space volume per parent energy shell:
        # g1/(2pi^2) int p E dE =g1/(6pi^2) [p_hi^3-p_lo^3].
        self.v1=g1/(6*math.pi**2)*(
            (ehi**2-M1_T*M1_T)**1.5-(elo**2-M1_T*M1_T)**1.5)
        self.pair_parent_index=np.add.outer(np.arange(nq),np.arange(nq))
        area=allowed_area_from_zero(self.edges,M1_T)
        valid=4*self.e[:,None]*self.e[None,:]>=M1_T*M1_T
        valid&=self.v1[self.pair_parent_index]>0
        self.K=(g1*M1_T/(2*math.pi**2))*area*valid
        self.active_pairs=int(np.count_nonzero(self.K))
        assert self.active_pairs>0

    def unpack(self,y):
        n=self.nq
        return y[:n],y[n:2*n],y[2*n:]

    def pair_flux(self,y):
        b,f,g=self.unpack(y)
        parent=g[self.pair_parent_index]
        return self.K*(parent*(1-f[None,:])*(1+b[:,None])-
                       (1-parent)*f[None,:]*b[:,None])

    def rhs(self,t,y):
        n=self.nq
        J=self.pair_flux(y)
        db=J.sum(axis=1)/self.vphi
        df=J.sum(axis=0)/self.vchi
        source=np.bincount(self.pair_parent_index.ravel(),weights=J.ravel(),
                           minlength=len(self.parent_E))
        dg=np.zeros(len(self.parent_E))
        good=self.v1>0
        dg[good]=-source[good]/self.v1[good]
        return np.r_[db,df,dg]

    def equilibrium(self,mu_phi=-.3,mu_chi=.1):
        if mu_phi>0:raise ValueError("Massless Bose mu cannot exceed zero")
        mu1=mu_phi+mu_chi
        return np.r_[1/np.expm1(self.e-mu_phi),
                      expit(mu_chi-self.e),expit(mu1-self.parent_E)]

    def integrals(self,y):
        b,f,g=self.unpack(y)
        Nphi=float(np.dot(self.vphi,b))
        Nchi=float(np.dot(self.vchi,f))
        N1=float(np.dot(self.v1,g))
        Ephi=float(np.dot(self.vphi*self.e,b))
        Echi=float(np.dot(self.vchi*self.e,f))
        E1=float(np.dot(self.v1*self.parent_E,g))
        return {
            "n_phi":Nphi,"n_chi2":Nchi,"n_chi1":N1,
            "charge_n_chi1_plus_chi2":N1+Nchi,
            "charge_n_chi1_plus_phi":N1+Nphi,
            "total_discrete_energy":E1+Echi+Ephi,
            "E_chi1":E1,"E_chi2":Echi,"E_phi":Ephi}

    def entropy_rate(self,y):
        b,f,g=self.unpack(y)
        if np.any(b<=0) or np.any(f<=0) or np.any(f>=1) or np.any(g<=0) or np.any(g>=1):
            raise ValueError("Quantum entropy diagnostic needs an interior state")
        d=self.rhs(0,y)
        n=self.nq
        return float(np.dot(self.vphi*d[:n],np.log(b/(1+b)))+
                     np.dot(self.vchi*d[n:2*n],np.log(f/(1-f)))+
                     np.dot(self.v1*d[2*n:],np.log(g/(1-g))))

    def entropy_production(self,y):
        b,f,g=self.unpack(y)
        if np.any(b<=0) or np.any(f<=0) or np.any(f>=1) or np.any(g<=0) or np.any(g>=1):
            raise ValueError("Interior only")
        a=g[self.pair_parent_index]*(1-f[None,:])*(1+b[:,None])
        v=(1-g[self.pair_parent_index])*f[None,:]*b[:,None]
        mask=self.K>0
        return float(np.sum((self.K*(a-v))[mask]*np.log(a[mask]/v[mask])))

    def checks(self):
        y=self.equilibrium(-.3,.1)
        d=self.rhs(0,y)
        equilibrium_error=float(np.max(np.abs(d*np.r_[self.vphi,self.vchi,self.v1])))
        init=np.r_[.1*np.exp(-self.e),.1*np.exp(-self.e),expit(-self.parent_E)]
        q=self.rhs(0,init)
        n=self.nq
        invariant1=float(np.dot(self.v1,q[2*n:])+np.dot(self.vchi,q[n:2*n]))
        invariant2=float(np.dot(self.v1,q[2*n:])+np.dot(self.vphi,q[:n]))
        energy=float(np.dot(self.v1*self.parent_E,q[2*n:])+
                     np.dot(self.vchi*self.e,q[n:2*n])+
                     np.dot(self.vphi*self.e,q[:n]))
        dH=self.entropy_rate(init)
        sigma=self.entropy_production(init)
        answer={"thermal_equilibrium_residual":equilibrium_error,
                "conserved_charge1_derivative":invariant1,
                "conserved_charge2_derivative":invariant2,
                "total_energy_derivative":energy,
                "negative_entropy_derivative":dH,
                "entropy_production":sigma,
                "dH_plus_sigma":dH+sigma}
        assert equilibrium_error<1e-11,answer
        assert max(abs(invariant1),abs(invariant2),abs(energy))<1e-10,answer
        assert dH<0 and sigma>0 and abs(dH+sigma)<1e-9,answer
        return answer

    def evolve(self,tend=.03):
        n=self.nq
        y0=np.r_[np.zeros(2*n),expit(-self.parent_E)]
        y0[2*n:][self.v1==0]=0
        initial=self.integrals(y0)
        sol=solve_ivp(self.rhs,(0,tend),y0,method="BDF",
                      rtol=1e-6,atol=1e-12,max_step=.0003)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1]
        b,f,g=self.unpack(y)
        final=self.integrals(y)
        residuals={key:float(final[key]-initial[key]) for key in (
             "charge_n_chi1_plus_chi2","charge_n_chi1_plus_phi",
             "total_discrete_energy")}
        assert max(map(abs,residuals.values()))<1e-8,residuals
        assert min(y)>=-1e-8 and max(f)<=1.000001 and max(g)<=1.000001
        return {"time_Gamma":tend,"initial":initial,"final":final,
                "conservation_residuals":residuals,
                "min_occupation":float(min(y)),
                "max_fermion_occupation_chi1":float(max(g)),
                "max_fermion_occupation_chi2":float(max(f)),
                "max_boson_occupation_grid_dependent":float(max(b)),
                "nfev":int(sol.nfev)}


def run(quick=False):
    configs=((30,.08),) if quick else ((30,.08),(50,.08),(80,.08),
                                       (100,.06),(150,.04))
    rows=[]
    for nq,delta in configs:
        k=ClosedResonantDecay(nq,delta)
        rows.append({"grid":{"nq":nq,"delta":delta,
                    "parent_energy_state_count":2*nq-1},
                    "checks":k.checks(),"evolution":k.evolve()})
    return {
      "status":"FLAT CLOSED THREE-SPECIES STRUCTURAL LATTICE TEST, NOT realistic IR production",
      "new_feature":"Parent chi1 dynamically depleted, not externally fixed",
      "daughter_energy_grid":"uniform; q_i=(i+.5)delta, parent E_{i+j}=(i+j+1)delta",
      "exact_discrete_properties":"two independent charges, total energy, Pauli/Bose positivity, detailed balance, H theorem",
      "limitations":"Finite uniform grid misses low-q burst, approximate cell rates and parent shell averages; no FRW, no rethermalizing SM bath, no physical freezeout, no CLASS prediction",
      "cases":rows}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/closed_three_species_lattice.json")
    opts=p.parse_args()
    result=run(opts.quick)
    out=Path(opts.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
