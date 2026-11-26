"""Exact momentum-differential VACUUM 1<->2 Bose-daughter source.

Thermal fermion parent at T=1, M/T=m, Gamma_rest=gamma, g1
helicity count; daughter phi/chi2 massless and initially EMPTY.

 dR/dq_phi = g1*m*gamma/(2*pi^2)
            log(1+exp(-q_phi-m^2/(4*q_phi))).

This peaks exactly at q=m/2 (NOT q=0); the deep infrared source
is exponentially suppressed before Bose stimulation develops.
This is an independent shape benchmark for external BEST C_phi(q):
  (g_phi/(2*pi²)) * q² * df_phi(q)/dt = dR/dq_phi, g_phi=1.

It is NOT an upper bound for nonlinear stimulated cold f(q), nor a
model of a condensate zero mode or cosmological freeze-out.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

from onshell_collision_2d import vacuum_reference


def vacuum_phi_source_differential(q, m_over_T=.1, g1=4,
                                   rest_width=1.):
    q=np.asarray(q,dtype=float)
    if np.any(~np.isfinite(q)) or np.any(q<0) or m_over_T<=0 or g1<=0 or rest_width<=0:
        raise ValueError("Valid nonnegative q, positive mass/dof/rate required")
    c=m_over_T*m_over_T/4
    exponent=np.divide(c,q,out=np.full_like(q,np.inf),where=q>0)
    spectrum=g1*m_over_T*rest_width/(2*np.pi**2) * np.logaddexp(
        0.,-q-exponent)
    if spectrum.ndim==0:return float(spectrum)
    return spectrum


def analytical_peak(m_over_T=.1,g1=4,rest_width=1.):
    q=m_over_T/2.
    return {"q_peak":q,
            "source_max_dR_dq":vacuum_phi_source_differential(q,m_over_T,g1,rest_width),
            "proof":"q+c/q >= 2sqrt(c)=m, equality at q=m/2; log(1+exp(-x)) monotone in -x",
            "kinematics":"M1/T=m, two massless daughters, initial thermal FD mother, empty daughters"}


def compare_best_phi_grid(q,phi_collision_rate,
                          m_over_T=.1,g1=4,rest_width=1.,gphi=1):
    momentum=np.asarray(q,dtype=float)
    source_rate=np.asarray(phi_collision_rate,dtype=float)
    if (momentum.ndim!=1 or len(momentum)<10 or
            source_rate.shape!=momentum.shape or
            np.any(~np.isfinite(momentum)) or
            np.any(~np.isfinite(source_rate)) or
            np.any(momentum<=0) or np.any(np.diff(momentum)<=0)
            or gphi<=0):
        raise ValueError("Valid strictly increasing positive BEST momentum grid and finite rates required")
    reference=vacuum_phi_source_differential(momentum,m_over_T,g1,rest_width)
    observed=gphi/(2*np.pi**2)*momentum**2*source_rate
    maximum=analytical_peak(m_over_T,g1,rest_width)["source_max_dR_dq"]
    selected=reference>maximum*.001
    if np.count_nonzero(selected)<3:
        raise ValueError("BEST grid does not cover enough of the physical vacuum source")
    residual=(observed[selected]-reference[selected])/maximum
    relative=(observed[selected]-reference[selected])/reference[selected]
    return {
        "q":momentum.tolist(),
        "BEST_dR_dq":observed.tolist(),
        "analytic_dR_dq":reference.tolist(),
        "fitted_source_q_peak":float(momentum[np.argmax(observed)]),
        "analytic_source_q_peak":m_over_T/2,
        "reference_relative_threshold_for_shape_comparison":.001,
        "robust_points_count":int(np.count_nonzero(selected)),
        "relative_peak_normalized_rms_difference":float(np.sqrt(np.mean(residual**2))),
        "maximum_absolute_relative_difference_within_support":float(np.max(abs(relative))),
        "max_absolute_source_error_as_fraction_of_analytic_peak":float(np.max(abs(residual))),
        "not_a_nonlinear_freezeout_shape_test":True,
        "warning":"Single initial Monte Carlo pass; finite q-grid and fluctuations. Vacuum source maximum q=m/2 is not a later Bose-stimulated peak or zero mode."
    }


def run():
    peak=analytical_peak()
    total=quad(lambda q:vacuum_phi_source_differential(q),
               0,np.inf,epsabs=1e-12,epsrel=1e-11)[0]
    analytic=vacuum_reference()
    if abs(total/analytic-1)>1e-10:
        raise AssertionError("Differential vacuum source fails number sum rule")
    samples=np.array([0.,1e-5,1e-4,.0005,.001,.005,.02,.05,.1,.5,1.,5.])
    return {
        "scope":"Exact initial differential source benchmark for EMPTY daughters only",
        "M_over_T":.1,"g1":4,"gamma":1.,
        "peak":peak,
        "number_rate_integrated_from_shape":float(total),
        "independent_1D_full_number_rate":float(analytic),
        "source_samples":[{"q":float(q),"dR_dq":float(v)}
                          for q,v in zip(samples,vacuum_phi_source_differential(samples))],
        "independent_BEST_executed":False,
        "warning":"Vacuum source may peak at q=M/2 even when nonlinear stimulated later solution produces vastly colder modes; do not extrapolate to CLASS."
    }


if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--output",default="outputs/vacuum_differential_spectrum.json")
    args=pa.parse_args()
    res=run()
    target=Path(args.output);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(res,indent=2),encoding="utf-8")
    print(json.dumps(res,indent=2))
