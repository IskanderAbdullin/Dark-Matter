"""On-shell 2D daughter-energy collision quadrature for chi1 <-> chi2 + phi.

ASSUMPTIONS: flat space, T=1, thermally FIXED Dirac mother with g1=4,
mass M1/T=0.1, massless daughters with g_phi=1, g_chi2=4, Gamma_rest=1.
No expansion, parent depletion or zero-mode condensate. The energy-pair
grid approximates the collision integral; it is NOT a cosmological model.

Change variables exactly:
 E1 = q_phi + q_chi,  p1 = sqrt(E1**2-M1**2);
 dmu/2 = dq_phi/p1, and 4*q_phi*q_chi>=M1**2.
 Thus (g1/2pi^2) p1^2 dp1 (M1/E1) (dmu/2)
 = [g1*M1/(2pi^2)] dq_phi dq_chi.
This eliminates mother and angular quadrature, sampling directly the
on-shell two-dimensional daughter-energy domain.

Cell areas of the hyperbolic domain are integrated ANALYTICALLY. The
positive collision kernel uses representative mean shell energies
(e_phi_i,e_chi_j), and drops partially allowed cells when their
representative pair violates the on-shell inequality. That introduces
a measured, vanishing-with-resolution threshold-cell discretization error.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.sparse import csr_matrix
from scipy.special import expit


def vacuum_reference(M1_T=.1,g1=4):
    """Continuum empty-daughter source, exactly integrated over q_chi.

    R = g1*M/(2pi^2) * integral_0^inf
        log(1+exp[-q-M^2/(4q)]) dq.
    """
    c=M1_T*M1_T/4
    def f(q):
        if q==0.:
            return 0.
        return float(np.logaddexp(0.,-q-c/q))
    a=quad(f,0.,1.,epsabs=1e-13,epsrel=1e-13,limit=250)[0]
    b=quad(f,1.,np.inf,epsabs=1e-13,epsrel=1e-13,limit=250)[0]
    return g1*M1_T/(2*math.pi**2)*(a+b)


def allowed_rectangle_areas(edges, M1_T=.1):
    """Exact area of rectangles intersected by q_phi*q_chi>=M^2/4.

    For a rectangular bin, integrate max(0,y_hi-max(y_lo,c/x)) dx:
    partial area = y_hi dx - c ln(x_b/x_a) for c/y_hi<x<c/y_lo,
    plus the fully allowed part x>=c/y_lo.
    """
    edges=np.asarray(edges,dtype=float)
    if np.any(edges<=0) or not np.all(np.diff(edges)>0):
        raise ValueError("strictly increasing positive edges required")
    c=M1_T*M1_T/4.
    xl=edges[:-1,None]
    xh=edges[1:,None]
    yl=edges[None,:-1]
    yh=edges[None,1:]
    xa=np.maximum(xl,c/yh)
    xb=np.minimum(xh,c/yl)
    part=np.where(xb>xa,yh*(xb-xa)-c*np.log(xb/xa),0.)
    xf=np.maximum(xl,c/yl)
    full=np.maximum(xh-xf,0.)*(yh-yl)
    area=np.maximum(part+full,0.)
    assert np.all(np.isfinite(area))
    assert np.max(area-(xh-xl)*(yh-yl))<1e-11
    return area


class OnShellCollision2D:
    def __init__(self,nq=150,qmin=1e-5,qmax=32.,M1_T=.1,g1=4,gchi=4):
        if nq<5 or qmin<=0 or qmax<=qmin or M1_T<=0:
            raise ValueError("Invalid grid or mass parameters")
        self.nq=nq;self.M1_T=M1_T;self.g1=g1;self.gchi=gchi
        self.edges=np.geomspace(qmin,qmax,nq+1)
        self.e=.75*(self.edges[1:]**4-self.edges[:-1]**4)/(self.edges[1:]**3-self.edges[:-1]**3)
        vol=(self.edges[1:]**3-self.edges[:-1]**3)/(6*math.pi**2)
        self.vol_phi=vol
        self.vol_chi=gchi*vol
        area=allowed_rectangle_areas(self.edges,M1_T)
        E=self.e[:,None]+self.e[None,:]
        kinetic_allowed=4*self.e[:,None]*self.e[None,:]>=M1_T*M1_T
        prefactor=g1*M1_T/(2*math.pi**2)
        full_kernel=prefactor*area
        thermal_a=expit(-E)
        raw_vac=float(np.sum(full_kernel*thermal_a))
        discarded=float(np.sum(full_kernel*thermal_a*(~kinetic_allowed)))
        # Every nonzero pair (e_i,e_j) is itself an on-shell 1<->2
        # reaction with E1=e_i+e_j, p1=sqrt(E1^2-M1^2), |mu|<=1.
        K=full_kernel*kinetic_allowed
        self.forward=csr_matrix(K*thermal_a)
        self.reverse=csr_matrix(K*(1-thermal_a))
        self.vacuum_rate=float(self.forward.sum())
        self.excluded_centroid_fraction=discarded/raw_vac if raw_vac else 0.
        self.vacuum_exact=vacuum_reference(M1_T,g1)
        self.vacuum_relative_error=self.vacuum_rate/self.vacuum_exact-1.
        self.pairs=int(np.count_nonzero(K))
        self.min_active_phi_q=float(self.e[np.any(K>0,axis=1)][0])
        self.max_allowed_pair_mass_deficit=float(np.max(
            np.maximum(M1_T*M1_T-4*self.e[:,None]*self.e[None,:],0.)*(K>0)))

    def rhs(self,t,y):
        b,f=y[:self.nq],y[self.nq:]
        db=(1+b)*(self.forward@(1-f))-b*(self.reverse@f)
        df=(1-f)*(self.forward.T@(1+b))-f*(self.reverse.T@b)
        return np.r_[db/self.vol_phi,df/self.vol_chi]

    def equilibrium(self,mu_phi=-.3):
        if mu_phi>0:raise ValueError("Bose chemical potential must be <=0")
        return np.r_[1/np.expm1(self.e-mu_phi),expit(-self.e-mu_phi)]

    def free_energy_rate(self,y):
        b,f=y[:self.nq],y[self.nq:]
        if np.any(b<=0) or np.any(f<=0) or np.any(f>=1):
            raise ValueError("Thermodynamic check requires interior state")
        d=self.rhs(0,y)
        return float(np.dot(self.vol_phi*d[:self.nq],self.e+np.log(b/(1+b)))+
                     np.dot(self.vol_chi*d[self.nq:],self.e+np.log(f/(1-f))))

    def entropy_production(self,y):
        b,f=y[:self.nq],y[self.nq:]
        if np.any(b<=0) or np.any(f<=0) or np.any(f>=1):
            raise ValueError("Thermodynamic check requires interior state")
        A=self.forward.tocoo()
        B=self.reverse.tocoo()
        assert np.array_equal(A.row,B.row) and np.array_equal(A.col,B.col)
        gain=A.data*(1+b[A.row])*(1-f[A.col])
        loss=B.data*b[B.row]*f[B.col]
        return float(np.dot(gain-loss,np.log(gain/loss)))

    def moment_budget(self,y):
        """Audit paired daughter number and on-shell parent reservoir energy.

        The parent is FIXED, not dynamically evolved. This is an explicit
        required energy transfer from the external thermal reservoir.
        """
        n=self.nq
        b,f=y[:n],y[n:]
        d=self.rhs(0,y)
        a=self.forward.tocoo()
        r=self.reverse.tocoo()
        assert np.array_equal(a.row,r.row) and np.array_equal(a.col,r.col)
        pair_rate=a.data*(1+b[a.row])*(1-f[a.col])-r.data*b[r.row]*f[r.col]
        pair_numbers=float(np.sum(pair_rate))
        pair_parent_energy=float(np.dot(pair_rate,self.e[a.row]+self.e[a.col]))
        nphi=float(np.dot(self.vol_phi,d[:n]))
        nchi=float(np.dot(self.vol_chi,d[n:]))
        Ephi=float(np.dot(self.e*self.vol_phi,d[:n]))
        Echi=float(np.dot(self.e*self.vol_chi,d[n:]))
        return {
            "reaction_number_rate":pair_numbers,
            "phi_number_rate":nphi,
            "chi2_number_rate":nchi,
            "daughter_number_difference_rate":nphi-nchi,
            "daughter_total_energy_rate":Ephi+Echi,
            "on_shell_parent_reservoir_energy_loss":-pair_parent_energy,
            "on_shell_energy_ledger_residual":Ephi+Echi-pair_parent_energy,
            "caveat":"Parent thermal distribution is fixed; only its required energy exchange is diagnosed."
        }

    def checks(self):
        d=self.rhs(0,self.equilibrium(-.3))
        equilibrium_error=float(max(np.max(abs(d[:self.nq]*self.vol_phi)),
                                    np.max(abs(d[self.nq:]*self.vol_chi))))
        y=np.r_[.02*np.exp(-self.e),.1*np.exp(-self.e)]
        dd=self.rhs(0,y)
        dF=self.free_energy_rate(y)
        sigma=self.entropy_production(y)
        z=np.zeros(2*self.nq)
        fone=z.copy();fone[self.nq:]=1.
        result={
          "equilibrium_max_bin_rate":equilibrium_error,
          "daughter_number_rate_mismatch":float(np.dot(dd[:self.nq],self.vol_phi)-
                                                np.dot(dd[self.nq:],self.vol_chi)),
          "free_energy_rate":dF,"entropy_production":sigma,"dF_plus_sigma":dF+sigma,
          "vacuum_boundary_valid":bool(np.min(self.rhs(0,z))>=-1e-12),
          "fermion_ceiling_valid":bool(np.max(self.rhs(0,fone)[self.nq:])<=1e-12)}
        assert equilibrium_error<1e-11,result
        assert abs(result["daughter_number_rate_mismatch"])<1e-10,result
        assert sigma>0 and dF<0 and abs(dF+sigma)<1e-9,result
        assert result["vacuum_boundary_valid"] and result["fermion_ceiling_valid"],result
        return result

    def evolve(self,tend=.03,rtol=1e-6,atol=1e-12,max_step=.0003):
        sol=solve_ivp(self.rhs,(0,tend),np.zeros(2*self.nq),method="BDF",
                      rtol=rtol,atol=atol,max_step=max_step)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1];b=y[:self.nq];f=y[self.nq:]
        mass=self.vol_phi*b
        nb=float(mass.sum())
        nchi=float(np.dot(self.vol_chi,f))
        cdf=np.cumsum(mass)/nb
        result={
          "time_Gamma":tend,
          "nphi_over_T3":nb,"nchi_over_T3":nchi,
          "daughter_difference_over_T3":nb-nchi,
          "mean_phi_q":float(np.dot(mass,self.e)/nb),
          "q10":float(np.interp(.1,cdf,self.e)),
          "q50":float(np.interp(.5,cdf,self.e)),
          "q90":float(np.interp(.9,cdf,self.e)),
          "fraction_phi_q_below_0p1":float(mass[self.e<.1].sum()/nb),
          "max_boson_occupation_grid_dependent":float(max(b)),
          "max_fermion_occupation":float(max(f)),
          "min_occupation":float(min(y)),
          "nfev":int(sol.nfev)}
        assert abs(nb-nchi)<2e-7,result
        assert result["min_occupation"]>=-2e-8,result
        assert result["max_fermion_occupation"]<=1.000001,result
        return result


def run(quick=False):
    grids=(50,100) if quick else (50,80,100,150,180,220,300)
    cases=[]
    for n in grids:
        k=OnShellCollision2D(n)
        cases.append({"nq":n,"vacuum_rate":k.vacuum_rate,
          "vacuum_relative_error":k.vacuum_relative_error,
          "exact_continuum_vacuum_rate":k.vacuum_exact,
          "excluded_invalid_centroid_fraction":k.excluded_centroid_fraction,
          "max_allowed_pair_mass_deficit":k.max_allowed_pair_mass_deficit,
          "active_pair_count":k.pairs,
          "checks":k.checks(),"evolved":k.evolve()})
    return {
      "description":"On-shell two-dimensional massless daughter 1<->2 quantum gain-loss",
      "baseline":"ea6ec461997600d5f00bf244ea808607d40ef8dd",
      "physical_validity":"Collision at representative daughter energies is on-shell; finite-grid area integration and centre threshold mask introduce controlled numerical error. Mother bath held fixed.",
      "continuum_vacuum_source_exact_integral":"g1 M1/(2pi²) int_0^infty dq log(1+exp[-q-M1²/(4q)])",
      "does_not_prove_FRW_frozen_spectrum":True,
      "scenarios":cases}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--quick",action="store_true")
    p.add_argument("--output",default="outputs/on_shell_2d_kinetics.json")
    opts=p.parse_args()
    ans=run(opts.quick)
    out=Path(opts.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(ans,indent=2),encoding="utf-8")
    print(json.dumps(ans,indent=2))
