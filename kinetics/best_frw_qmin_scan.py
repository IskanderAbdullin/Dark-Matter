"""Independent BEST radiation-FRW qmin-volume scan, not a BEC inference.

At a fixed REAL time resolution (128 steps), q-grid (96 bins),
Vegas neval and physical model, vary only the COMOVING smallest q.
The full quantum 3-species + massive zero-mode code is called from
best_radiation_frw_pilot.py (not forked or modified).

Compare measured BEST zero-mode N0(qmin) against the *independently
integrated* regular massive finite-q Bose particle number
  N_IR(qmin)=gphi/(2pi²) ∫_0^qmin q² f(q,t)dq
using physically redshifted original FD mother initial condition.

If N0 ∝ qmin³ and N0/qmin³ stabilizes, it represents ordinary
smooth IR phase-volume, not macroscopic delta-function condensate.
A true BEC test requires late-time qmin and grid/seed/time convergence.

No execution of external BEST without explicit --execute.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from best_radiation_frw_pilot import rd_protocol,execute_best
from best_massive_daughter_reference import massive_protocol
from radiation_collisionless_parent_bose import evaluate


def qmin_protocol(qmin=1e-4,grid_n=96,neval=6000,
                  steps=128,dt=3.125e-7,t0=1e-4):
    p=rd_protocol(t0=t0,steps=steps,dtGamma=dt,nq=grid_n,neval=neval)
    p["kinetic_source"]=massive_protocol(
        mass=.1,mphi=.001,mchi=0.,g1=4,rest_width=1.,
        qmin=qmin,qmax=30.,grid_n=grid_n,vegas_neval=neval)
    p["purpose"]="Independent massive φ FRW zero-mode qmin phase-volume convergence at fixed end a"
    p["qmin_scan_is_not_macroscopic_BEC"]=True
    p["nonzero_BEC_condensate_inferred"]=False
    p["external_BEST_executed"]=False
    return p


def execute_case(best_root,output,p):
    original=execute_best(best_root,output,p)
    if original is None:return None
    qmin=p["kinetic_source"]["qmin"]
    t0=p["rd_t0"]
    t_end=t0+p["end_tGamma"]
    estimate=evaluate(t0=t0,end_t=t_end,qmin=qmin)
    last=original["time_series"][-1]
    original["independent_regular_phase_volume_reference"]=estimate
    original["zero_mode_measured_to_analytic_regular_cell_ratio"]=(
        last["phi_massive_zero_mode_N0"]/
        estimate["below_qmin_comoving_number"])
    original["not_automatically_a_condensate"]=True
    original["independent_BEST_qmin_scan_requires_multiple_same_qcut_grids"]=True
    Path(output).write_text(json.dumps(original,indent=2),encoding="utf-8")
    print(json.dumps({
        "qmin":qmin,
        "measured_zero_mode_N0":last["phi_massive_zero_mode_N0"],
        "analytic_regular_IR_number":estimate["below_qmin_comoving_number"],
        "zero_mode_ratio":original["zero_mode_measured_to_analytic_regular_cell_ratio"],
        "measured_phi_f_qmin":last["phi_occupation_at_smallest_positive_q"],
        "analytic_phi_f_qmin":estimate["fphi_at_qmin"],
        "not_convergence_or_BEC_claim":True
    },indent=2))
    return original


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--execute",action="store_true")
    p.add_argument("--best-root",default="")
    p.add_argument("--output",default="outputs/best_frw_qmin.json")
    p.add_argument("--qmin",type=float,default=1e-4)
    p.add_argument("--grid-n",type=int,default=96)
    p.add_argument("--neval",type=int,default=6000)
    p.add_argument("--steps",type=int,default=128)
    p.add_argument("--dt",type=float,default=3.125e-7)
    opts=p.parse_args()
    protocol=qmin_protocol(opts.qmin,opts.grid_n,opts.neval,opts.steps,opts.dt)
    if opts.execute:
        if not opts.best_root:raise ValueError("--best-root required")
        execute_case(opts.best_root,opts.output,protocol)
    else:
        print(json.dumps(protocol,indent=2))


if __name__=="__main__":
    main()
