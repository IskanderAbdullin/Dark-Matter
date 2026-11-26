"""First combined physically normalized five-population Dirac/anti real-Bose
Yukawa quantum WEAK operator: both decays and BOTH pair channels.

Species order: chi1,anti_chi1,chi2,anti_chi2,real_phi.
- χ1 -> χ2+φ / antiχ1 -> antiχ2+φ use the identical boosted on-shell
  physical parent/angle events previously derived in the coupled 3-field
  Yukawa model, splitting original g1=4 rate into g_particle=2 each.
- χ1+antiχ1 -> 2φ via virtual χ2, and χ2+antiχ2 -> 2φ via virtual χ1,
  use independently derived positive physical 5D lab COM scattering events.
All four channels share each channel's OWN positive forward/reverse W,
exact energy shell and the SAME physical Bose occupations at actual E.

This is an instantaneous continuum quantum collision weak source. It does
NOT yet provide a finite-volume 5-field implicit strong solver, in-medium
thermal dispersion, SM portal or self-consistent FRW relic.
"""
import json
import math
import numpy as np
from scipy.special import expit

from kinetics.phi_bose_coupled_bernstein_p2_finite import Coupled
from kinetics.phi_bose_pair_physical_onshell_quantum import PhysicalPairEvents
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI

SPECIES=("chi1","anti_chi1","chi2","anti_chi2","phi")
CHARGES=np.array([1.,-1.,1.,-1.,0.])

def _eta(energies,temperature,mu,curv=0.):
    if temperature<=0:raise ValueError("Positive physical temperatures required")
    return (energies-mu)/temperature+curv*(1-np.cos(energies))

def _occ(eta,kind):
    if kind=="fermion":return expit(-eta)
    if np.min(eta)<=0:raise ValueError("Nonpositive Bose inverse fugacity")
    return 1/np.expm1(eta)

class FullDiracYukawaWeak:
    """A physical shared-event four-reaction network over five populations."""
    def __init__(self,P=2.,Ppair=None,nphi=12,nrad=8,nangle=12,
                 pair_nrad=8,pair_nrelative=8,pair_ncm=32,pair_naz=4,
                 pair_radial_breaks=None):
        self.P=P
        self.Ppair=P if Ppair is None else Ppair
        self.decay=Coupled(P=P,Nphi=nphi,nrad=nrad,nangle=nangle)
        kw=dict(P=self.Ppair,nrad=pair_nrad,nrelative=pair_nrelative,
                ncm=pair_ncm,naz=pair_naz,radial_breaks=pair_radial_breaks)
        self.pair1=PhysicalPairEvents(
            **kw,m_in=M1,m_virtual=M2)
        self.pair2=PhysicalPairEvents(
            **kw,m_in=M2,m_virtual=M1)
        self.event_counts={
            "chi1_decay":len(self.decay.W),
            "anti_chi1_decay":len(self.decay.W),
            "chi1_pair":len(self.pair1.W),
            "chi2_pair":len(self.pair2.W)}

    def source(self,T=(.9,.9,1.,1.,1.),mu=(0.,0.,0.,0.,0.),
               curvature=(0.,0.,0.,0.,0.),
               include_decay=True,include_pairs=True):
        if len(T)!=5 or len(mu)!=5 or len(curvature)!=5:
            raise ValueError("Exactly five Dirac/anti and real-scalar states")
        T=np.asarray(T,float)
        mu=np.asarray(mu,float)
        curv=np.asarray(curvature,float)
        N=np.zeros(5)
        Q=np.zeros(5)
        H=0.
        rates={}
        events=[]
        def shared_event(name,species,E,kind,weights):
            nonlocal H,N,Q
            etas=[_eta(e,T[i],mu[i],curv[i]) for i,e in zip(species,E)]
            f=[_occ(e,k) for e,k in zip(etas,kind)]
            if len(species)==3:
                # species=(parent, daughter, phi), event sign parent-1
                i,j,k=species
                F=f[0]*(1-f[1])*(1+f[2])
                R=(1-f[0])*f[1]*f[2]
                stoich=np.array([-1.,1.,1.])
            else:
                # species=(Dirac particle, antidirac, phi3, phi4)
                i,j,k,l=species
                F=f[0]*f[1]*(1+f[2])*(1+f[3])
                R=(1-f[0])*(1-f[1])*f[2]*f[3]
                stoich=np.array([-1.,-1.,1.,1.])
            J=F-R
            signed=J*weights
            aff=sum(s*eta for s,eta in zip(stoich,etas))
            h=float(np.dot(signed,aff))
            # This check may permit tiny machine-level cancellation
            if min(float(np.min(J*aff)),0.) < -1e-11:
                raise ArithmeticError("Physical quantum microscopic H sign")
            H+=h
            rate=float(np.sum(signed))
            for i,s,e in zip(species,stoich,E):
                N[i]+=s*rate
                Q[i]+=s*float(np.dot(signed,e))
            rates[name]=dict(event_rate_GeV4=rate,quantum_H_GeV4=h,
                             min_micro_JA=float(np.min(J*aff)),
                             min_weight=float(np.min(weights)),
                             event_count=len(weights))
        E=self.decay.E
        if include_decay:
            for name,i,j in (("chi1_decay",0,2),
                             ("anti_chi1_decay",1,3)):
                shared_event(name,(i,j,4),(E[0],E[1],E[2]),
                             ("fermion","fermion","boson"),
                             self.decay.W/2)
        if include_pairs:
            for name,model,i,j in (
                ("chi1_pair",self.pair1,0,1),
                ("chi2_pair",self.pair2,2,3)):
                shared_event(name,(i,j,4,4),
                             (model.e_chi,model.e_anti,model.e_phi3,model.e_phi4),
                             ("fermion","fermion","boson","boson"),
                             model.W)
        return dict(dN=N.tolist(),dQ=Q.tolist(),H_GeV4=H,
                    vector_U1_defect_GeV4=float(CHARGES@N),
                    total_energy_defect_GeV5=float(np.sum(Q)),
                    event_channels=rates,
                    class_ready=False)

    def audit(self):
        out={}
        out["equilibrium"]=self.source(T=(1.,)*5,
                                        mu=(.05,-.05,.05,-.05,0.))
        out["nonthermal"]=self.source(
            T=(.9,1.05,1.,1.1,.98),
            mu=(0.,0.,0.,0.,0.),
            curvature=(.02,.01,.015,-.01,.005))
        out["restricted_decay_mu_not_full_equilibrium"]=self.source(
            T=(1.,)*5,mu=(-.15,.15,.05,-.05,-.2))
        out["model"]=dict(P=self.P,Ppair=self.Ppair,
                          event_counts=self.event_counts,
                          included_charges=CHARGES.tolist())
        out["class_ready"]=False
        return out

if __name__=="__main__":
    m=FullDiracYukawaWeak()
    print(json.dumps(m.audit(),indent=2,allow_nan=False))
