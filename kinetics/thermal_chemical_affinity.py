"""Exact thermal-occupancy 1<->2 net source at arbitrary chemical affinity.

Massive FD mother chi1 -> massless FD chi2 + massless BE phi, with
same temperature T=1 but arbitrary (constant) chemical potentials.
All daughter momenta are physical on-shell if 4*q_phi*q_chi >= m**2.

For q=q_phi, y=q_chi, y0=m²/(4q),
  A=e^(q-mu1), B=e^(-mu2), t0=e^(-y0)
  int_{y0}^infty dy fFD(y+q-mu1)*(1-fFD(y-mu2))
     = B/(A-B) * [log1p(t0/B) - log1p(t0/A)].

Multiply by (1+bphi) and by
  1-exp(mu2+mu_phi-mu1) = -expm1(-chemical_affinity).
This EXACTLY integrates net quantum gain-minus-loss in y.
Use a stable A≈B limit and the zero source at q=0.

The sign of the integrated net source is therefore exactly the sign
of mu1-mu2-mu_phi, for all q>0 in this common-T thermal snapshot.
This is an instantaneous nonlinear quantum collision benchmark,
not a solved time-dependent FRW/condensate/dark-matter relic.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from best_external_vacuum_reference import benchmark_protocol


def thermal_1to2_spectrum(q, mass=.1, g1=4, rest_width=1.,
                         mu1=-.15, mu2=.1, mu_phi=-.3):
    q=np.asarray(q,dtype=float)
    if (mass<=0 or g1<=0 or rest_width<=0 or mu_phi>=0
            or np.any(~np.isfinite(q)) or np.any(q<0)
            or not all(np.isfinite(x) for x in (mu1,mu2,mu_phi))):
        raise ValueError("Positive mass, g1, width; mu_phi<0; finite q and chemical potentials required")
    affinity=mu1-mu2-mu_phi
    x=q.reshape(-1)
    out=np.zeros_like(x)
    indices=np.flatnonzero((x>0)&(x<100))
    if len(indices):
        p=x[indices]
        threshold=mass*mass/(4*p)
        t=np.exp(-threshold)
        A=np.exp(p-mu1)
        B=np.exp(-mu2)
        diff=A-B
        I=np.empty_like(p)
        near=np.abs(diff)<1e-7*np.maximum(A,B)
        I[near]=t[near]/(A[near]+t[near])
        far=~near
        I[far]=B/diff[far]*(np.log1p(t[far]/B)-np.log1p(t[far]/A[far]))
        Bphi=1/np.expm1(p-mu_phi)
        forward=(1+Bphi)*I
        out[indices]=g1*mass*rest_width/(2*np.pi**2)*forward*(-np.expm1(-affinity))
    if out.size==1 and q.ndim==0:
        return float(out[0])
    return out.reshape(q.shape)


def integrated_thermal_source(mass=.1,g1=4,rest_width=1.,
                              mu1=-.15,mu2=.1,mu_phi=-.3):
    """At common T=1, source n/T^3/dt; all q integrated."""
    keyword=dict(mass=mass,g1=g1,rest_width=rest_width,
                 mu1=mu1,mu2=mu2,mu_phi=mu_phi)
    total=quad(lambda q:thermal_1to2_spectrum(q,**keyword),0,1.,
               epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    total+=quad(lambda q:thermal_1to2_spectrum(q,**keyword),1,np.inf,
                epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    energy_phi=quad(lambda q:q*thermal_1to2_spectrum(q,**keyword),
                    0,1,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    energy_phi+=quad(lambda q:q*thermal_1to2_spectrum(q,**keyword),
                     1,np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    return {
        "mass_over_T":mass,"mu1":mu1,"mu2":mu2,"mu_phi":mu_phi,
        "chemical_affinity":mu1-mu2-mu_phi,
        "n_phi_net_rate":float(total),
        "n_chi2_net_rate":float(total),
        "n_chi1_net_rate":float(-total),
        "phi_net_physical_energy_rate":float(energy_phi),
        "net_rate_sign_matches_affinity":bool(
            total==0 if mu1-mu2-mu_phi==0
            else total*(mu1-mu2-mu_phi)>0),
        "not_a_time_evolution":True,
    }


def independent_2d_integrand(q,y,mass=.1,mu1=-.15,mu2=.1,mu_phi=-.3):
    """Direct Bose/Pauli forward-backward difference, no y integration."""
    if q<=0 or y<=0 or 4*q*y<mass*mass:
        raise ValueError("Physical on-shell daughter energies required")
    a=float(expit(-(q+y-mu1)))
    f=float(expit(-(y-mu2)))
    b=1/np.expm1(q-mu_phi)
    return float(a*(1-f)*(1+b)-(1-a)*f*b)


def independent_physical_energy_reference(mass=.1,g1=4,rest_width=1.,
                                          mu1=-.15,mu2=.1,mu_phi=-.3):
    """Independent original 2D gain-loss quadrature for all three energy rates.

    Unlike the closed-form daughter y-integration this evaluates the
    ACTUAL Bose/Pauli terms inside two nested adaptive quadratures.
    Upper q=50 has exponentially negligible tail for |chemical mu|<=5.
    It is a fixed-temperature source ledger, not FRW pressure work.
    """
    if max(abs(mu1),abs(mu2),abs(mu_phi))>5:
        raise ValueError("Finite 50T energy quadrature requires |mu|<=5")
    _=thermal_1to2_spectrum(.1,mass,g1,rest_width,mu1,mu2,mu_phi)
    pref=g1*mass*rest_width/(2*np.pi**2)
    def nested(q,kind):
        if q<=0:return 0.
        ymin=mass*mass/(4*q)
        def integrand(y):
            net=independent_2d_integrand(q,y,mass,mu1,mu2,mu_phi)
            E=q if kind=="phi" else (y if kind=="chi2" else q+y)
            return E*net
        return quad(integrand,ymin,np.inf,epsabs=1e-10,
                    epsrel=1e-9,limit=200)[0]
    energies={}
    for species,sign in (("phi",1),("chi2",1),("chi1",-1)):
        flux=quad(lambda q:nested(q,species),0,1,epsabs=1e-10,
                  epsrel=1e-9,limit=200)[0]
        flux+=quad(lambda q:nested(q,species),1,50,epsabs=1e-10,
                   epsrel=1e-9,limit=200)[0]
        energies[species]=sign*pref*flux
    analytic_phi=integrated_thermal_source(
        mass,g1,rest_width,mu1,mu2,mu_phi)["phi_net_physical_energy_rate"]
    residual=sum(energies.values())
    if abs(residual)>1e-8 or abs(energies["phi"]-analytic_phi)>1e-8:
        raise ArithmeticError("Direct 2D energy or analytical phi source is inconsistent")
    return {
        "physical_energy_derivatives":energies,
        "energy_conservation_residual":residual,
        "analytic_phi_energy_rate":analytic_phi,
        "analytic_vs_independent_2d_phi_energy_error":energies["phi"]-analytic_phi,
        "two_dimensional_original_quantum_factors_used":True,
        "q_upper_cutoff_for_this_reference":50.,
        "strictly_instantaneous_1to2_source":True
    }


def run():
    cases=[integrated_thermal_source(mu1=value)
           for value in (-.25,-.2,-.15)]
    vacuum=integrated_thermal_source(mu1=0,mu2=-30.,mu_phi=-30.)
    target=benchmark_protocol()["independently_derived_vacuum_source_n_over_T3_per_time"]
    assert abs(vacuum["n_phi_net_rate"]/target-1)<1e-9
    return {
      "model":"exact Bose-enhanced and Pauli-blocked thermal 1<->2 source at fixed common T",
      "mass_over_T":.1,"g1":4,"rest_width":1,
      "cases":cases,
      "empty_daughter_chemical_limit":vacuum,
      "independent_empty_daughter_number_source":target,
      "vacuum_relative_error":vacuum["n_phi_net_rate"]/target-1,
      "external_BEST_executed":False,
      "warning":"Instantaneous chemical-affinity source only, no full nonlinear time evolution or cold DM spectrum."
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="outputs/thermal_chemical_affinity.json")
    args=parser.parse_args()
    data=run()
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(data,indent=2),encoding="utf-8")
    print(json.dumps(data,indent=2))
