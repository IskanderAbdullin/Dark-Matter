"""Convergent low-q Bose source cell versus the exact q=0 occupation seed.

For an initially empty, finite-mass phi produced in a thermal 1->2
decay, the exact differentiated rate is dR/dq = (g_phi/2pi²) q² df/dt.
The occupation at the exact point q=0 may be HUGE while the number
added to the physical phase-volume cell 0<q<qcut is negligible.

Compute the ACTUAL number flux into this cell using an independent
adaptive quadrature of the finite-mass on-shell differential source:
 R_cell=integral_0^qcut dR/dq dq,
 f_cell_dot = R_cell/[gphi*qcut³/(6pi²)].

The special Bose zero mode/condensate n0 is NOT the entire finite-q
cell, and this initial vacuum source is not evidence of nonlinear BEC.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from massive_daughter_decay_source import (
    massive_daughter_vacuum_spectrum,
    zero_momentum_bose_occupation_source,
)
from onshell_collision_2d import vacuum_reference


def initial_ir_cell(q_cut=.0001,M=.1,mphi=.001,mchi=0.,
                    g1=4,gphi=1,rest_width=1.):
    if q_cut<=0 or gphi<=0 or not math.isfinite(q_cut):
        raise ValueError("Positive finite q cutoff and Bose degeneracy required")
    kw={"M":M,"mphi":mphi,"mchi":mchi,"g1":g1,"rest_width":rest_width}
    result=quad(lambda x:massive_daughter_vacuum_spectrum(x,**kw),
                0,q_cut,epsabs=1e-17,epsrel=1e-10,limit=350)[0]
    if result<0:
        raise ArithmeticError("Negative initial vacuum Bose number flux")
    phase_vol=gphi*q_cut**3/(6*math.pi**2)
    mean_source=result/phase_vol
    f0=zero_momentum_bose_occupation_source(
        M=M,mphi=mphi,mchi=mchi,g1=g1,gphi=gphi,
        rest_width=rest_width)
    full_source=vacuum_reference(M,g1)*rest_width
    return {
        "M_over_T":M,"mphi_over_T":mphi,"mchi_over_T":mchi,
        "q_cut":q_cut,"gphi":gphi,"Gamma_rest":rest_width,
        "qcell_particle_phase_volume":phase_vol,
        "exact_physical_number_flux_below_qcut":float(result),
        "exact_cell_phase_volume_averaged_occupation_rate":float(mean_source),
        "exact_qzero_occupation_rate":float(f0),
        "cell_average_to_qzero_source_ratio":
            float(mean_source/f0) if f0>0 else None,
        "fraction_all_initial_phi_decays_below_qcut":
            float(result/full_source),
        "is_a_condensate_n0_source":False,
        "physical_meaning":"Finite positive-q phase volume. q=0 occupation source alone is NOT a delta-function condensate particle number."
    }


def run():
    rows=[initial_ir_cell(mphi=m,q_cut=cut)
          for m in (.0005,.001,.002,.005)
          for cut in (1e-5,5e-5,1e-4)]
    return {
        "scope":"Exact initial vacuum Bose source below finite q cutoff for massive phi",
        "cases":rows,
        "important_point":"Extremely large finite q=0 occupation sources can coexist with a <1e-7 fraction of all produced particles below qcut=1e-4.",
        "independent_BEST_run_required_for_massive_phi":True,
        "zero_mode_or_FRW_classification":"NOT VALIDATED: requires explicit n0 kinetics, self-consistent time evolution and energy/entropy budgets"
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/massive_phi_ir_cell.json")
    o=p.parse_args()
    d=run()
    file=Path(o.output);file.parent.mkdir(parents=True,exist_ok=True)
    file.write_text(json.dumps(d,indent=2),encoding="utf-8")
    print(json.dumps(d,indent=2))
