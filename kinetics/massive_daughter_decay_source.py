"""EXACT finite-mass two-body thermal mother decay spectrum, no Bose backreaction.

A massive thermal fermion chi1 (M=M1/T), constant rest width Gamma,
decays isotropically to phi (mphi/T) and chi2 (mchi/T). Initial
daughters are empty. The rest daughter energy and momentum are

 Ephi*=(M²+mphi²-mchi²)/(2M),
 p*=sqrt([M²-(mphi+mchi)²][M²-(mphi-mchi)²])/(2M).

A daughter of laboratory momentum q has E=sqrt(q²+mphi²).
From energy-momentum conservation the mother energy is restricted to
  Emin=M*(Ephi*²+q²)/(E*E+q*p*)
  Emax=M*(E*E+q*p*)/mphi², if mphi>0 (else infinity).
The exact vacuum source at T=1 is

 dR/dq = g1*Gamma*M²/(4*pi²*p*) * q/E *
         [log(1+exp(-Emin))-log(1+exp(-Emax))].

Use a cancellation-safe difference of log1p values near q=0.

At FIXED physical Gamma and an OPEN decay, integrated daughter number
equals the usual parent thermal time-dilated decay rate INDEPENDENT
of daughter masses. Energy fractions are Ephi*/M and Echi*/M, but
the IR mode-dependent differential occupation is NOT continuous under
setting the physical boson mass to zero. For any mphi>0:

  (d fphi(0)/dt) = (g1/gphi)*Gamma*(M/mphi)^3 *
                   fFD(M*Ephi*/mphi).

This conditional seed may be huge even if the massless approximation
looks good for the total number. The zero-mode phase-volume integration
and subsequent Bose stimulation are NOT implemented here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from onshell_collision_2d import vacuum_reference
from frw_reservoir_budget import thermal_mother_number


def rest_kinematics(M=.1,mphi=0.,mchi=0.):
    if not all(math.isfinite(v) for v in (M,mphi,mchi)) or M<=0 or mphi<0 or mchi<0:
        raise ValueError("Positive mother mass and nonnegative finite daughter masses required")
    if mphi+mchi>=M:
        raise ValueError("Threshold closed; fixed nonzero rest width is not physical")
    estar=(M*M+mphi*mphi-mchi*mchi)/(2*M)
    estar_chi=(M*M+mchi*mchi-mphi*mphi)/(2*M)
    pstar=math.sqrt((M*M-(mphi+mchi)**2)*(M*M-(mphi-mchi)**2))/(2*M)
    return {
        "M1_over_T":M,"mphi_over_T":mphi,"mchi_over_T":mchi,
        "rest_phi_energy":estar,"rest_chi_energy":estar_chi,
        "rest_daughter_momentum":pstar,
        "phi_energy_fraction":estar/M,
        "chi_energy_fraction":estar_chi/M,
    }


def massive_daughter_vacuum_spectrum(q,M=.1,mphi=0.,mchi=0.,
                                     g1=4,rest_width=1.):
    k=rest_kinematics(M,mphi,mchi)
    if g1<=0 or rest_width<=0:
        raise ValueError("Positive mother degeneracy and decay width required")
    momentum=np.asarray(q,dtype=float)
    if not np.all(np.isfinite(momentum)) or np.any(momentum<0):
        raise ValueError("Finite nonnegative daughter momenta required")
    flat=momentum.reshape(-1)
    values=np.zeros_like(flat)
    positive=flat>0
    if np.any(positive):
        Q=flat[positive]
        E=np.hypot(Q,mphi)
        estar=k["rest_phi_energy"]
        pstar=k["rest_daughter_momentum"]
        x=E*estar+Q*pstar
        lower=M*(estar**2+Q**2)/x
        if mphi>0:
            upper=M*x/(mphi*mphi)
            gap=2*M*Q*pstar/(mphi*mphi)
            # Stable even if upper-lower is smaller than float rounding.
            gain=np.log1p(np.exp(-lower)*(-np.expm1(-gap))/
                           (1+np.exp(-upper)))
        else:
            gain=np.log1p(np.exp(-lower))
        pref=g1*rest_width*M*M/(4*math.pi**2*pstar)
        values[positive]=pref*(Q/E)*gain
    if momentum.ndim==0:
        return float(values[0])
    return values.reshape(momentum.shape)


def zero_momentum_bose_occupation_source(M=.1,mphi=.001,mchi=0.,
                                         g1=4,gphi=1,rest_width=1.):
    """Exact finite q->0 occupation SOURCE, not a condensate density."""
    kin=rest_kinematics(M,mphi,mchi)
    if g1<=0 or gphi<=0 or rest_width<=0:
        raise ValueError("Positive degeneracy and width required")
    if mphi==0:
        return 0.
    eparent=M*kin["rest_phi_energy"]/mphi
    log_value=(math.log(g1/gphi)+math.log(rest_width)+
               3*math.log(M/mphi)-float(np.logaddexp(0.,eparent)))
    return float(math.exp(log_value)) if log_value>-745 else 0.


def integrated_source_and_energy(M=.1,mphi=0.,mchi=0.,g1=4,
                                 rest_width=1.):
    kin=rest_kinematics(M,mphi,mchi)
    def spectrum(q):
        return massive_daughter_vacuum_spectrum(
            q,M,mphi,mchi,g1,rest_width)
    N=quad(spectrum,0,1,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    N+=quad(spectrum,1,np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    phi_energy=quad(lambda q:math.hypot(q,mphi)*spectrum(q),
                    0,1,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    phi_energy+=quad(lambda q:math.hypot(q,mphi)*spectrum(q),
                     1,np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    exact_N=vacuum_reference(M,g1)*rest_width
    exact_total_energy=M*rest_width*thermal_mother_number(M,g1)
    expected_phi_energy=kin["phi_energy_fraction"]*exact_total_energy
    return {
        "kinematics":kin,
        "total_phi_number_rate":N,
        "fixed_gamma_exact_number_rate":exact_N,
        "number_integral_relative_error":N/exact_N-1,
        "direct_phi_physical_energy_rate":phi_energy,
        "rest_frame_exact_phi_physical_energy_rate":expected_phi_energy,
        "phi_physical_energy_relative_error":phi_energy/expected_phi_energy-1,
        "fixed_gamma_parent_physical_energy_loss":exact_total_energy,
        "exact_chi_energy_fraction":kin["chi_energy_fraction"],
        "finite_mass_q0_Bose_occupation_source":
            zero_momentum_bose_occupation_source(M,mphi,mchi,g1,1,rest_width),
        "fixed_width_assumption":"Gamma_rest held fixed even when daughter masses vary; this is NOT the fixed-matrix-element threshold dependence",
        "no_Bose_stimulation_or_zero_mode_evolution":True,
    }


def run():
    ratios=(0.,.001,.005,.01,.02,.05,.1,.2)
    cases=[{"mphi_over_M":r,**integrated_source_and_energy(mphi=.1*r)}
           for r in ratios]
    assert all(abs(x["number_integral_relative_error"])<1e-9 for x in cases)
    assert all(abs(x["phi_physical_energy_relative_error"])<1e-9 for x in cases)
    return {
        "description":"Exact finite-mass isotropic two-body vacuum daughter spectrum and q=0 occupation source",
        "M1_over_T":.1,"g1":4,"gphi":1,"rest_width":1.,
        "cases":cases,
        "important_limit":"Nonzero massive phi has potentially very large finite q=0 occupation SOURCE; strictly massless phi source tends to zero at q=0. Small mphi does not mean uniformly small error in IR.",
        "research_scope":"Instantaneous vacuum-source-only, fixed Gamma and mother thermal FD. No nonlinear condensate, evolving population, SM plasma, FRW or CLASS."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/finite_mass_daughter_source.json")
    opts=p.parse_args()
    result=run()
    target=Path(opts.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))
