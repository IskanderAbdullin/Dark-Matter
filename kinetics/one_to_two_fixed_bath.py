"""Conservative fixed-bath 1<->2 Bose stimulation benchmark.

Implements the local statistical gain/loss for chi1 <-> chi2 + phi as in
Wen Yin, JHEP05(2023)180, eqs. (2.2)-(2.4),(2.9)-(2.11), section 3.
Only Minkowski space, fixed T=1, fixed thermal mother chi1 bath, massless
daughter chi2 and massless boson phi. M1/T=.1 and eta=1. Gamma_rest=1.
Does NOT solve FRW expansion, mother-bath depletion, realistic dark-sector
particle Lagrangian or the Bose-star/galactic kinetic equations.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit,roots_legendre

class FlatThermalDecay:
    def __init__(self,nq=150,np1=180,nmu=96,M1_T=.1):
        self.nq=nq
        self.M1_T=M1_T
        self.edges=np.geomspace(1e-5,32.,nq+1)
        # One spin state phase-space volume q^2 dq/(2 pi^2), q measured in T.
        volume=(self.edges[1:]**3-self.edges[:-1]**3)/(6*math.pi**2)
        self.phi_volume=volume
        self.chi_volume=4*volume
        self.qcentre=np.sqrt(self.edges[1:]*self.edges[:-1])
        self.avg_q=(3./4)*(self.edges[1:]**4-self.edges[:-1]**4)/(self.edges[1:]**3-self.edges[:-1]**3)
        # chi1 is a thermal Dirac fermion with g1=4. Two-body decay
        # yields Ephi=(E1+p1*cos(theta))/2, Echi2=E1-Ephi.
        pn,pw=roots_legendre(np1)
        p=13*(pn+1)
        wp=13*pw
        mu,wm=roots_legendre(nmu)
        P=p[:,None]
        E=np.sqrt(P*P+M1_T*M1_T)
        qphi=(E+P*mu[None,:])/2
        qchi=E-qphi
        self.qphi=qphi.ravel()
        self.qchi=qchi.ravel()
        self.parentE=np.broadcast_to(E,qphi.shape).ravel()
        self.parent_fermi=np.broadcast_to(expit(-E),qphi.shape).ravel()
        # mother density: g1/(2pi^2) p1^2 dp; decay gamma_lab=M1/E1;
        # isotropic daughter rest-frame dmu/2; common event rate for
        # both daughters guarantees exact *number* conservation.
        self.rateprefactor=(4/(2*math.pi**2)*np.broadcast_to(
            P**2*wp[:,None]*(M1_T/E)*(wm[None,:]/2),qphi.shape)).ravel()
        self.ip,self.wp=self._interp_weights(self.qphi)
        self.ic,self.wc=self._interp_weights(self.qchi)

    def _interp_weights(self,q):
        logcentre=np.log(self.qcentre)
        lo=np.clip(np.searchsorted(self.qcentre,q)-1,0,self.nq-2)
        t=(np.log(q)-logcentre[lo])/(logcentre[lo+1]-logcentre[lo])
        return lo,np.clip(t,0.,1.)

    def event_gain_loss(self,y,stimulated=True):
        fp=y[:self.nq]
        fc=y[self.nq:]
        phi=fp[self.ip]*(1-self.wp)+fp[self.ip+1]*self.wp
        chi=fc[self.ic]*(1-self.wc)+fc[self.ic+1]*self.wc
        a=self.parent_fermi
        # Both parent and chi2 are fermions. phi is a boson.
        S=a*(1-chi)*((1+phi) if stimulated else 1.)-(1-a)*chi*phi
        return self.rateprefactor*S

    def rhs(self,t,y,stimulated=True):
        rates=self.event_gain_loss(y,stimulated)
        dotphi=(np.bincount(self.ip,weights=rates*(1-self.wp),minlength=self.nq)
            +np.bincount(self.ip+1,weights=rates*self.wp,minlength=self.nq))/self.phi_volume
        dotchi=(np.bincount(self.ic,weights=rates*(1-self.wc),minlength=self.nq)
            +np.bincount(self.ic+1,weights=rates*self.wc,minlength=self.nq))/self.chi_volume
        return np.r_[dotphi,dotchi]

    def detailed_balance_error(self):
        # Exact equilibrium occupations at EVENT momenta, not interpolated
        # bin-centre values. Energy relation is exact for every event.
        fphi=1/np.expm1(self.qphi)
        fchi=expit(-self.qchi)
        a=self.parent_fermi
        S=a*(1-fchi)*(1+fphi)-(1-a)*fchi*fphi
        return max(abs(S)),max(abs(self.qphi+self.qchi-self.parentE))

    def snapshot(self,y,time,stimulated):
        fphi=y[:self.nq];fchi=y[self.nq:]
        dphi=fphi*self.phi_volume
        dchi=fchi*self.chi_volume
        nphi=float(dphi.sum())
        nchi=float(dchi.sum())
        return {"time_Gamma_rest":time,"stimulated":bool(stimulated),
          "nphi_over_T3":nphi,"nchi_over_T3":nchi,
          "daughter_number_mismatch":float((nphi-nchi)/max(nphi,1e-200)),
          "phi_fraction_with_q_below_0p1":float(dphi[self.qcentre<.1].sum()/max(nphi,1e-200)),
          "phi_mean_q":float(np.dot(dphi,self.avg_q)/max(nphi,1e-200)),
          "phi_rms_q":float(np.sqrt(np.dot(dphi,self.avg_q**2)/max(nphi,1e-200))),
          "phi_max_occupation_not_pointwise_converged":float(max(fphi)),
          "chi2_max_fermi_occupation":float(max(fchi)),
          "min_occupation":float(min(min(fphi),min(fchi)))}

    def evolve(self,stimulated,tend=.03):
        ts=[.0001,.001,.003,.01,tend]
        ts=sorted(set(t for t in ts if t<=tend))
        y0=np.zeros(2*self.nq)
        sol=solve_ivp(lambda t,y:self.rhs(t,y,stimulated),
          [0,tend],y0,t_eval=ts,method="RK45",
          max_step=.0006,rtol=1e-5,atol=2e-10)
        if not sol.success:
            raise RuntimeError(sol.message)
        records=[self.snapshot(y,t,stimulated) for y,t in zip(sol.y.T,sol.t)]
        for r in records:
            assert abs(r["daughter_number_mismatch"])<1e-8,r
            assert r["min_occupation"]>-1e-8,r
            assert r["chi2_max_fermi_occupation"]<1.00001,r
        return {"nfev":sol.nfev,"history":records,"final":records[-1]}

def run(quick=False):
    grids=([(100,140,72)] if quick else
       [(100,140,72),(150,180,96),(180,240,120)])
    results=[]
    for nq,np1,nmu in grids:
        kernel=FlatThermalDecay(nq,np1,nmu)
        balance,energy=kernel.detailed_balance_error()
        assert balance<1e-10 and energy<1e-12,(balance,energy)
        bose=kernel.evolve(True)
        control=kernel.evolve(False)
        b=bose["final"];n=control["final"]
        assert b["nphi_over_T3"]>100*n["nphi_over_T3"]
        assert b["phi_fraction_with_q_below_0p1"]>.99
        assert n["phi_fraction_with_q_below_0p1"]<.15
        results.append({
          "grid":{"nq":nq,"mother_quadrature":np1,"decay_angle_quadrature":nmu,
                  "event_cells":np1*nmu},
          "thermal_event_detailed_balance_max_abs":float(balance),
          "event_energy_conservation_max_abs":float(energy),
          "with_Bose_stimulation":bose,
          "without_stimulation_diagnostic_only":control,
          "number_amplification_final":b["nphi_over_T3"]/n["nphi_over_T3"]})
        print("grid",nq,np1,nmu,"Bose density",b["nphi_over_T3"],
           "no-stimulation density",n["nphi_over_T3"],"cold fraction",b["phi_fraction_with_q_below_0p1"],
           flush=True)
    if len(results)>1:
        nums=np.array([x["with_Bose_stimulation"]["final"]["nphi_over_T3"] for x in results])
        mom=np.array([x["with_Bose_stimulation"]["final"]["phi_mean_q"] for x in results])
        densityspread=float(nums.max()/nums.min()-1)
        momentumspread=float(mom.max()/mom.min()-1)
        assert densityspread<.03,("integrated density did not converge",densityspread)
        shape_is_resolved=momentumspread<.05
    else:
        densityspread=None;momentumspread=None;shape_is_resolved=False
    return {
      "source":"Wen Yin JHEP05(2023)180 collision factor and 1<->2 kinematics",
      "exact_assumptions":"M1/T=.1,M2=0,mphi~0 during production; chi1 Dirac thermally FIXED, chi2 Dirac kinetic, phi real Bose kinetic; flat space, Gamma_rest=1",
      "raw_reaction":"S=f1*(1-f2)*(1+fphi)-(1-f1)*f2*fphi",
      "occupation_interpolation":"log-momentum cloud-in-cell; conservative daughter numbers, thermal event-level detailed balance",
      "grids":results,
      "relative_integrated_phi_number_spread":densityspread,
      "relative_mean_phi_momentum_spread":momentumspread,
      "momentum_shape_converged_to_5pct":shape_is_resolved,
      "warning":"Integrated number converges but momentum shape and maximum mode occupancy do not. Do NOT use this spectrum as a physical CLASS initial PSD yet. FRW redshift, source depletion, source rethermalization, collision between species, inverse cooling and number at fixed cosmological entropy remain untreated."
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/kinetic_1to2_bose_burst.json")
    args=parser.parse_args()
    out=run(args.quick)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print("saved",path)
