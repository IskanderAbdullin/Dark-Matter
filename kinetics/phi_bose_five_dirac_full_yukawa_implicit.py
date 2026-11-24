"""Conditional fully implicit five-population Yukawa FD/BE physical moment prototype.

Microscopic network and real on-shell continuum event weights:
  chi1 -> chi2 + phi, anti_chi1 -> anti_chi2 + phi,
  chi1 + anti_chi1 -> phi + phi, chi2 + anti_chi2 -> phi + phi.
Both 1↔2 and 2↔2 processes use ONE positive common quantum F/R event
measure at each event. Nonnegative degree-2 Bernstein basis reproduces
physical E EXACTLY, with g*p²dp/(2π²) per species physical DOS capacity.
Separate particles and antiparticles have spin g=2 each, real phi g=1.
Only net Dirac U1, total on-shell energy, H theorem are fundamental
linear invariants of this four-channel five-field graph.

A positive solved fully implicit moment root obeys a conditional
quantum relative-entropy finite-step inequality. This is a NUMERICAL
DARK-SECTOR vacuum collision prototype: no in-medium χ/φ dispersion,
scalar quartic, SM portal, Friedmann, zero-mode, freezeout, Ω or CLASS.
No unconditional existence/solver convergence or pointwise eta-tangent
positivity from Bernstein-P2 is asserted.
"""
import math
import numpy as np
from scipy.special import expit
from scipy.sparse import csr_matrix,block_diag
from scipy.sparse.linalg import spsolve

from kinetics.phi_bose_coupled_bernstein_p2_finite import (
    build_species,bernstein
)
from kinetics.phi_bose_five_dirac_full_yukawa_weak import (
    FullDiracYukawaWeak,CHARGES,SPECIES
)
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI
from kinetics.phi_bose_pair_ir_parent_support import global_pair_scalar_support_floor

def asflat(x):return np.concatenate(x)
def unflat(m,x):return [np.asarray(x[a:b]) for a,b in
                         zip(m.offsets[:-1],m.offsets[1:])]

class EventChannel:
    """Each event leg matrix A reconstructs actual η(E), S reproduces E."""
    def __init__(self,name,indices,stoich,energies,weights,model):
        self.name=name
        self.indices=tuple(indices)
        self.stoich=np.asarray(stoich,float)
        self.W=np.asarray(weights,float)
        self.energy=tuple(np.asarray(e,float) for e in energies)
        n=len(self.W)
        if len(self.indices)!=len(self.stoich) or len(self.indices)!=len(self.energy):
            raise ValueError("Stoichiometry and microscopic event legs disagree")
        if not n or np.min(self.W)<=0 or not all(len(e)==n for e in self.energy):
            raise ValueError("Nonpositive physical Yukawa event table")
        self.kind=[model.species[i].kind for i in self.indices]
        self.A=[]
        for i,e in zip(self.indices,self.energy):
            s=model.species[i]
            c=np.searchsorted(s.edges,e,side="right")-1
            c=np.clip(c,0,s.ncells-1)
            L=s.edges[c];U=s.edges[c+1]
            x=(e-L)/(U-L)
            if np.min(x)<-4e-9 or np.max(x)>1+4e-9:
                raise ValueError(("Physical event out of physical species DOS grid",name,i,
                                  float(np.min(x)),float(np.max(x))))
            hat=bernstein(np.clip(x,0,1))
            rows=np.repeat(np.arange(n),3)
            cols=(model.offsets[i]+3*c[:,None]+np.arange(3)).ravel()
            self.A.append(csr_matrix((hat.ravel(),(rows,cols)),
                              shape=(n,model.nvars)))
        self.S=sum((sgn*A for sgn,A in zip(self.stoich,self.A)),
                   csr_matrix((n,model.nvars)))
        self.max_energy_defect=float(np.max(np.abs(
           sum(float(sgn)*e for sgn,e in zip(self.stoich,self.energy)))))
        if self.max_energy_defect>1e-12:
            raise ArithmeticError("Invalid four-leg physical on-shell energy")
        # The same exact sparse basis generates the correct physical
        # four-leg affinity and all number/energy ledgers.
        self.model=model

    def weights_and_derivatives(self,eta,with_jac=False):
        etas=[a@eta for a in self.A]
        occ=[]
        for e,k in zip(etas,self.kind):
            if k=="boson":
                if np.min(e)<=0:raise ValueError("Nonpositive Bose η")
                occ.append(1/np.expm1(e))
            else:occ.append(expit(-e))
        F=np.ones_like(self.W)
        R=np.ones_like(self.W)
        for f,kind,sgn in zip(occ,self.kind,self.stoich):
            if kind=="fermion":
                if sgn<0:F*=f;R*=1-f
                else:F*=1-f;R*=f
            else:
                if sgn!=1:raise ValueError("Only outgoing real Bose in listed Yukawa channels")
                F*=1+f;R*=f
        J=F-R
        affinity=self.S@eta
        H=float(np.dot(self.W*J,affinity))
        out=dict(J=J,F=F,R=R,H=H,
                 event_rate=float(np.dot(self.W,J)),
                 minimum_event_entropy=float(np.min(J*affinity)))
        if not with_jac:return out
        dJ=[]
        for f,kind,sgn in zip(occ,self.kind,self.stoich):
            if kind=="fermion":
                dj=(-F*(1-f)-R*f) if sgn<0 else (F*f+R*(1-f))
            else:
                dj=-F*f+R*(1+f)
            dJ.append(dj)
        out["dJ"]=dJ
        return out


class FullDiracPhysicalMoments:
    """Five-field physical DOS Bernstein-P2 strong finite-time model."""
    def __init__(self,P=2.,Ppair=None,N1=3,N2=3,Nphi=8,
                 nq=24,nrad=4,nangle=8,pair_nrad=4,
                 pair_nrelative=8,pair_ncm=32,pair_naz=4):
        self.weak=FullDiracYukawaWeak(
            P=P,Ppair=Ppair,nphi=Nphi,nrad=nrad,nangle=nangle,
            pair_nrad=pair_nrad,pair_nrelative=pair_nrelative,
            pair_ncm=pair_ncm,pair_naz=pair_naz)
        self.P=P
        self.Ppair=self.weak.Ppair
        e1hi=max(math.hypot(self.Ppair,M1),
                 self.weak.decay.species[0].edges[-1])
        e2hi=max(math.hypot(self.Ppair,M2),
                 self.weak.decay.species[1].edges[-1])
        # Pair COM boosted daughter energy may approach the sum of incoming
        # Echi+Eanti, not the old decay-only φ upper shell.
        phihi=2*math.hypot(self.Ppair,max(M1,M2))
        # The lighter Dirac-antipair channel reaches a much lower
        # physically allowed massive-scalar energy than 1↔2 decay.
        # Split at BOTH exact hard-P pair support floors, so active pair
        # source cannot leak through a broad first Bose entropy-dual cell
        # into truly inaccessible φ energies below the first allowed event.
        pair_minima=[global_pair_scalar_support_floor(self.Ppair,m)
                     for m in (M1,M2)]
        phi_edges=np.unique(np.r_[self.weak.decay.species[2].edges,
                                pair_minima,
                                max(phihi,1.05*self.weak.decay.species[2].edges[-1])])
        s1=build_species("chi1",M1,2,"fermion",np.linspace(M1,e1hi,N1+1),nq)
        sa1=build_species("anti_chi1",M1,2,"fermion",
                           np.linspace(M1,e1hi,N1+1),nq)
        s2=build_species("chi2",M2,2,"fermion",np.linspace(M2,e2hi,N2+1),nq)
        sa2=build_species("anti_chi2",M2,2,"fermion",
                           np.linspace(M2,e2hi,N2+1),nq)
        sp=build_species("real_phi",MPHI,1,"boson",phi_edges,nq)
        self.species=[s1,sa1,s2,sa2,sp]
        self.offsets=np.cumsum([0]+[s.nvars for s in self.species])
        self.nvars=int(self.offsets[-1])
        d=self.weak.decay
        self.channels=[
            EventChannel("chi1_decay",(0,2,4),(-1,+1,+1),
                         (d.E[0],d.E[1],d.E[2]),d.W/2,self),
            EventChannel("anti_chi1_decay",(1,3,4),(-1,+1,+1),
                         (d.E[0],d.E[1],d.E[2]),d.W/2,self)]
        for name,px,i,j in (
            ("chi1_pair",self.weak.pair1,0,1),
            ("chi2_pair",self.weak.pair2,2,3)):
            self.channels.append(EventChannel(
                name,(i,j,4,4),(-1,-1,+1,+1),
                (px.e_chi,px.e_anti,px.e_phi3,px.e_phi4),px.W,self))
        self.total_events=sum(len(c.W) for c in self.channels)
        self.charge_nodal=asflat([
            np.full(s.nvars,CHARGES[i]) for i,s in enumerate(self.species)])
        self.energy_nodal=asflat([s.nodal_energies for s in self.species])
        for ch in self.channels:
            if np.max(abs(ch.S@self.charge_nodal))>5e-12:
                raise ArithmeticError("Net real Dirac U1 broken by physical event")
            if np.max(abs(ch.S@self.energy_nodal))>5e-12:
                raise ArithmeticError("Physical on-shell energy broken by event interpolation")

    def initial_etas(self,T=(.9,1.05,1.,1.1,.98),
                     mu=(0.,0.,0.,0.,0.)):
        return [(s.nodal_energies-mu[i])/T[i]
                for i,s in enumerate(self.species)]

    def source(self,etas,with_jac=False):
        vec=asflat(etas)
        C=np.zeros(self.nvars)
        H=0.
        M=[]
        terms=[]
        for ch in self.channels:
            x=ch.weights_and_derivatives(vec,with_jac=with_jac)
            C+=np.asarray(ch.S.T@(ch.W*x["J"])).ravel()
            H+=x["H"]
            terms.append(dict(name=ch.name,event_rate=x["event_rate"],
                              H=x["H"],min_JA=x["minimum_event_entropy"]))
            if with_jac:
                dC=sum((ch.S.T@a.multiply((ch.W*dj)[:,None])
                        for a,dj in zip(ch.A,x["dJ"])),
                       csr_matrix((self.nvars,self.nvars)))
                M.append(dC)
        out=dict(C=C,H=H,channels=terms,
                 exact_event_U1=float(self.charge_nodal@C),
                 exact_event_energy=float(self.energy_nodal@C))
        if with_jac:
            out["dC"]=sum(M,csr_matrix((self.nvars,self.nvars)))
        return out

    def integrate_spectrum(self,etas):
        vals=[s.evaluate(x,with_entropy=True) for s,x in zip(self.species,etas)]
        N=np.asarray([sum(x["M"]) for x in vals])
        Q=np.asarray([s.nodal_energies@x["M"] for s,x in zip(self.species,vals)])
        S=sum(x["S"] for x in vals)
        return dict(N=N,Q=Q,S=S,U1=CHARGES@N)

    def implicit_residual(self,new,old,dt,with_jac=True):
        vec=asflat(unflat(self,new))
        prev=asflat(unflat(self,old))
        es=unflat(self,vec)
        os=unflat(self,prev)
        delta=asflat([s.moment_difference(e,o) for s,e,o in
                      zip(self.species,es,os)])
        source=self.source(es,with_jac=with_jac)
        R=delta-dt*source["C"]
        if not with_jac:return R,source
        G=block_diag([g for s,e in zip(self.species,es)
                      for g in s.evaluate(e,with_G=True)["G"]],
                     format="csr")
        Jac=-G-dt*source["dC"]
        return R,Jac,source

    def implicit_step(self,etas,dt,maxiter=24,rtol=2e-7):
        if dt<=0:raise ValueError("Positive dt required")
        old=asflat(etas)
        initial=self.integrate_spectrum(etas)
        src0=self.source(etas)
        flat_M=asflat([s.evaluate(e)["M"] for s,e in zip(self.species,etas)])
        scale=np.maximum.reduce([dt*np.abs(src0["C"]),
                    1e-14*np.abs(flat_M),np.full(self.nvars,1e-24)])
        if max(abs(dt*src0["C"])/scale)<1e-9:
            return dict(etas=etas,iterations=0,scaled_residual=0.,
                        entropy_increase=0.,entropy_certificate=0.,
                        energy_defect=0.,U1_defect=0.,physical_H=src0["H"],
                        min_Bose_eta=float(min(etas[4])),
                        stalled_at_precision=True,accepted_dt=dt)
        vec=old.copy()
        errors=[]
        for i in range(maxiter):
            R,J,info=self.implicit_residual(vec,old,dt)
            err=float(max(abs(R/scale)))
            errors.append(err)
            if err<rtol:break
            update=spsolve(J,-R)
            if not np.all(np.isfinite(update)):
                raise ArithmeticError("Non-finite 5-field sparse Newton update")
            alpha=1.
            accepted=False
            for _ in range(32):
                trial=vec+alpha*update
                if np.min(trial[self.offsets[4]:])<=0:
                    alpha/=2;continue
                try:
                    residual,_=self.implicit_residual(trial,old,dt,False)
                except (ArithmeticError,ValueError,OverflowError,FloatingPointError):
                    alpha/=2;continue
                next_err=float(max(abs(residual/scale)))
                if next_err<=err*(1-1e-4*alpha) or next_err<rtol:
                    vec=trial
                    accepted=True
                    break
                alpha/=2
            if not accepted:raise ArithmeticError(("5field Newton line search stalled",
                                                    dt,i,err))
        else:
            raise ArithmeticError(("5field Newton failed",dt,errors[-1]))
        updated=unflat(self,vec)
        end=self.integrate_spectrum(updated)
        s=self.source(updated)
        dS=end["S"]-initial["S"]
        energy=float(sum(end["Q"]-initial["Q"]))
        u1=float(end["U1"]-initial["U1"])
        return dict(etas=updated,iterations=i,scaled_residual=errors[-1],
                    entropy_increase=dS,entropy_certificate=dS-dt*s["H"],
                    energy_defect=energy,U1_defect=u1,physical_H=s["H"],
                    min_Bose_eta=float(min(updated[4])),
                    source_new=s,accepted_dt=dt,residual_history=errors)

    def adaptive_step(self,etas,dt,maxhalvings=12):
        failures=[]
        for _ in range(maxhalvings+1):
            try:
                out=self.implicit_step(etas,dt)
                scale=max(1e-25,dt*sum(abs(x["event_rate"]) for x in
                            self.source(etas)["channels"]))
                if out["entropy_increase"]<-1e-8*scale:
                    raise ArithmeticError("Negative entropy")
                if out["entropy_certificate"]<-1e-8*scale:
                    raise ArithmeticError("Negative quantum entropy gap")
                if max(abs(out["energy_defect"]),abs(out["U1_defect"]))>1e-8*scale:
                    raise ArithmeticError("Physical U1 or energy deficit")
                if out["min_Bose_eta"]<=0:
                    raise ArithmeticError("Bose positivity lost")
                out["rejections"]=failures
                return out
            except (ArithmeticError,ValueError,np.linalg.LinAlgError,OverflowError) as e:
                failures.append((dt,str(e)))
                dt/=2
        raise ArithmeticError(("5field all trial steps rejected",failures))

    def audit(self):
        test={}
        for tag,etas in (
            ("nonthermal",self.initial_etas()),
            ("thermal_full_U1",self.initial_etas(T=(1.,)*5,
                                              mu=(.05,-.05,.05,-.05,0.))),
            ("thermal_mu_phi_negative",self.initial_etas(
                 T=(1.,)*5,mu=(-.15,.15,.05,-.05,-.2)))):
            src=self.source(etas)
            test[tag]=dict(source_H=src["H"],event_energy=src["exact_event_energy"],
                           event_U1=src["exact_event_U1"],
                           channels=src["channels"])
            if tag=="nonthermal":
                a=self.adaptive_step(etas,1e5)
                test[tag]["finite_step"]={k:v for k,v in a.items()
                                         if k not in ("etas","source_new")}
        return dict(schema_version=1,scope="first fully coupled 5-population Yukawa physical-capacity quantum FD/BE implicit finite collision step; NOT SM/FRW/86eV relic",
                    P=self.P,Ppair=self.Ppair,nvars=self.nvars,
                    events=self.total_events,measurements=test,
                    class_ready=False)
if __name__=="__main__":
    import json
    result=FullDiracPhysicalMoments().audit()
    print(json.dumps(result,indent=2,default=lambda x:x.tolist() if hasattr(x,"tolist") else str(x)))
