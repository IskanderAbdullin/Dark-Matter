"""Finite-grid entropy-balanced 1<->2 SURROGATE, not a physical FRW collision solver.

Use the existing parent quadrature/decay angles, but evaluate mother thermal
occupation at SUM OF RECIPIENT SHELL ENERGIES, not its original physical E1.
This makes a finite lattice gain/loss exactly detailed-balanced at thermal
BE/FD occupations. The approximation error to physical parent E1 is measured.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from scipy.special import expit
from one_to_two_fixed_bath import FlatThermalDecay


class EntropyBalancedLattice:
    def __init__(self, nq=100, np1=140, nmu=72, qmin=1e-5, angular_rule="mu"):
        if not (0 < qmin <= 1e-5):
            raise ValueError("qmin outside supported (0,1e-5] interval")
        self.angular_rule=angular_rule
        if angular_rule=="mu":
            k=FlatThermalDecay(nq=nq,np1=np1,nmu=nmu)
        elif angular_rule=="logq":
            from logq_angular_quadrature import LogQDaughterKernel
            k=LogQDaughterKernel(nq=nq,np1=np1,nmu=nmu)
        else:
            raise ValueError("angular_rule must be mu or logq")
        self.nq, self.qmin = nq, qmin
        if qmin != 1e-5:
            k.edges=np.geomspace(qmin,32.,nq+1)
            k.phi_volume=(k.edges[1:]**3-k.edges[:-1]**3)/(6*math.pi**2)
            k.chi_volume=4*k.phi_volume
            k.qcentre=np.sqrt(k.edges[:-1]*k.edges[1:])
            k.avg_q=.75*(k.edges[1:]**4-k.edges[:-1]**4)/(k.edges[1:]**3-k.edges[:-1]**3)
            k.ip,k.wp=k._interp_weights(k.qphi)
            k.ic,k.wc=k._interp_weights(k.qchi)
        self.e=k.avg_q
        self.vol_b=k.phi_volume
        self.vol_f=k.chi_volume
        self.qcentre=k.qcentre
        self.min_quadrature_daughter_q=float(min(k.qphi.min(),k.qchi.min()))
        self.angular_normalization_error=getattr(k,"angular_normalization_max_error",None)
        rows=[];cols=[];fwd=[];rev=[]
        true_vac=float(np.dot(k.rateprefactor,k.parent_fermi))
        abs_vac_err=0.;abs_ener_err=0.;weighted_ener=0.
        for i in (0,1):
            ip=k.ip+i
            wp=k.wp if i else 1-k.wp
            for j in (0,1):
                ic=k.ic+j
                wc=k.wc if j else 1-k.wc
                weight=k.rateprefactor*wp*wc
                Eeff=self.e[ip]+self.e[ic]
                ae=expit(-Eeff)
                rows.append(ip);cols.append(ic)
                fwd.append(weight*ae);rev.append(weight*(1-ae))
                abs_vac_err+=float(np.dot(weight,np.abs(ae-k.parent_fermi)))
                abs_ener_err+=float(np.dot(weight,np.abs(Eeff-k.parentE)))
                weighted_ener+=float(np.dot(weight,k.parentE))
        ij=(np.concatenate(rows),np.concatenate(cols))
        self.forward=coo_matrix((np.concatenate(fwd),ij),shape=(nq,nq)).tocsr()
        self.backward=coo_matrix((np.concatenate(rev),ij),shape=(nq,nq)).tocsr()
        self.true_vacuum_rate=true_vac
        self.lattice_vacuum_rate=float(self.forward.sum())
        self.relative_vacuum_shift=self.lattice_vacuum_rate/true_vac-1.
        self.absolute_forward_weight_error=abs_vac_err/true_vac
        self.absolute_parent_energy_error=abs_ener_err/weighted_ener

    def rhs(self, t, y):
        b,f=y[:self.nq],y[self.nq:]
        db=(1+b)*(self.forward@(1-f))-b*(self.backward@f)
        df=(1-f)*(self.forward.T@(1+b))-f*(self.backward.T@b)
        return np.r_[db/self.vol_b,df/self.vol_f]

    def equilibrium(self, mu_b=-0.3):
        if mu_b>0:
            raise ValueError("mu_b cannot be positive for massless BE gas")
        return np.r_[1/np.expm1(self.e-mu_b),expit(-self.e-mu_b)]

    def free_energy_rate(self, y):
        b,f=y[:self.nq],y[self.nq:]
        if not (np.all(b>0) and np.all((f>0)&(f<1))):
            raise ValueError("Free-energy calculation requires interior state")
        d=self.rhs(0,y)
        return float(np.dot(self.vol_b*d[:self.nq],self.e+np.log(b/(1+b)))+
                     np.dot(self.vol_f*d[self.nq:],self.e+np.log(f/(1-f))))

    def entropy_production(self, y):
        b,f=y[:self.nq],y[self.nq:]
        if not (np.all(b>0) and np.all((f>0)&(f<1))):
            raise ValueError("Entropy calculation requires interior state")
        A=self.forward.tocoo()
        B=self.backward.tocoo()
        assert np.array_equal(A.row,B.row) and np.array_equal(A.col,B.col)
        gain=A.data*(1+b[A.row])*(1-f[A.col])
        loss=B.data*b[B.row]*f[B.col]
        return float(np.dot(gain-loss,np.log(gain/loss)))

    def checks(self):
        eq=self.equilibrium(-0.3)
        d=self.rhs(0,eq)
        residual=float(max(np.max(np.abs(d[:self.nq]*self.vol_b)),
                           np.max(np.abs(d[self.nq:]*self.vol_f))))
        x=np.r_[0.02*np.exp(-self.e),0.1*np.exp(-self.e)]
        dx=self.rhs(0,x)
        dF=self.free_energy_rate(x)
        sigma=self.entropy_production(x)
        z=np.zeros(2*self.nq)
        fone=z.copy();fone[self.nq:]=1
        return {
            "thermal_equilibrium_max_bin_rate":residual,
            "daughter_number_rate_mismatch":float(np.dot(dx[:self.nq],self.vol_b)-
                                                  np.dot(dx[self.nq:],self.vol_f)),
            "free_energy_derivative":dF,
            "entropy_production":sigma,
            "dF_plus_entropy_production":dF+sigma,
            "vacuum_inward":bool(np.min(self.rhs(0,z))>=-1e-12),
            "fermion_ceiling_inward":bool(np.max(self.rhs(0,fone)[self.nq:])<=1e-12),
        }

    def evolve(self, tend=0.03):
        sol=solve_ivp(self.rhs,(0,tend),np.zeros(2*self.nq),
                      method="BDF",rtol=1e-6,atol=1e-12,max_step=.0003)
        if not sol.success:
            raise RuntimeError(sol.message)
        y=sol.y[:,-1];b=y[:self.nq];f=y[self.nq:]
        nphi=float(np.dot(self.vol_b,b))
        nchi=float(np.dot(self.vol_f,f))
        result={
            "t_Gamma":tend, "nphi_over_T3":nphi,"nchi_over_T3":nchi,
            "daughter_number_difference_over_T3":nphi-nchi,
            "mean_phi_q":float(np.dot(self.vol_b*b,self.e)/nphi),
            "phi_fraction_below_q0p1":float((self.vol_b*b)[self.e<0.1].sum()/nphi),
            "max_boson_occupation_grid_dependent":float(np.max(b)),
            "max_fermion_occupation":float(np.max(f)),
            "min_occupation":float(np.min(y)),
            "nfev":int(sol.nfev)}
        assert abs(nphi-nchi)<2e-7,result
        assert result["min_occupation"]>=-2e-8,result
        assert result["max_fermion_occupation"]<=1.000001,result
        return result


def scenario(nq,np1,nmu,qmin=1e-5,integrate=True,angular_rule="mu"):
    model=EntropyBalancedLattice(nq,np1,nmu,qmin=qmin,angular_rule=angular_rule)
    checks=model.checks()
    assert checks["thermal_equilibrium_max_bin_rate"]<1e-11,checks
    assert abs(checks["daughter_number_rate_mismatch"])<1e-10,checks
    assert checks["free_energy_derivative"]<0 and checks["entropy_production"]>0,checks
    assert abs(checks["dF_plus_entropy_production"])<1e-9,checks
    assert checks["vacuum_inward"] and checks["fermion_ceiling_inward"],checks
    return {
        "grid":[nq,np1,nmu],"qmin":qmin,"angular_rule":angular_rule,
        "lowest_quadrature_daughter_q":model.min_quadrature_daughter_q,
        "angular_normalization_error":model.angular_normalization_error,
        "true_parent_vacuum_rate":model.true_vacuum_rate,
        "lattice_surrogate_vacuum_rate":model.lattice_vacuum_rate,
        "relative_total_vacuum_rate_shift":model.relative_vacuum_shift,
        "relative_absolute_forward_rate_projection_error":model.absolute_forward_weight_error,
        "relative_absolute_parent_energy_projection_error":model.absolute_parent_energy_error,
        "structural_checks":checks,
        "evolution":model.evolve() if integrate else None}


def run(quick=False,angular_rule="mu"):
    configs=[(50,64,48)] if quick else [
        (50,64,48),(100,140,72),(150,180,96),(180,240,120)]
    grid=[scenario(*dims,angular_rule=angular_rule) for dims in configs]
    ir=[] if quick else [scenario(150,180,96,qmin=q,angular_rule=angular_rule) for q in (1e-6,1e-7)]
    return {
        "model":"Entropy-balanced finite lattice with EFFECTIVE mother energy at recipient-pair centres",
        "angular_rule":angular_rule,
        "baseline":"ea6ec461997600d5f00bf244ea808607d40ef8dd",
        "physical_parent_energy_not_exact":True,
        "not_a_frw_or_class_initial_condition":True,
        "grid_study":grid,
        "infrared_cutoff_study":ir,
        "warning":"Mother E1 in Fermi factor is replaced by e_phi_bin+e_chi_bin. Discrete balance/entropy apply only to the surrogate. Gamma, expansion, source depletion, dynamic condensate and physical freeze-out are not modelled."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--angular-rule",choices=("mu","logq"),default="mu")
    parser.add_argument("--output",default="outputs/entropy_balanced_lattice.json")
    args=parser.parse_args()
    output=run(args.quick,args.angular_rule)
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(output,indent=2,ensure_ascii=False))
