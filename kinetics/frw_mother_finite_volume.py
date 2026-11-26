"""Finite-volume massive-mother FRW Liouville moments on fixed q=p/T bins.

Occupations f_i are the EXACT phase-volume-averaged initial FD state
(16-point Gauss quadrature) in each comoving q bin. They remain constant
without collisions. This guarantees number conservation at the
semidiscrete level without artificial mother rethermalization.

For each subsequent mass/temperature x=M1/T, compute rho/P/d rho/dx
by positive quadrature using the SAME bin f_i. The exact algebraic
identity x d(rho/T^4)/dx = rho/T^4-3P/T^4 holds node by node.

WARNING: even perfect number conservation and FRW continuity do NOT
prove the BIN approximation converged for energy/pressure. Test against
independent continuous FD moments and vary the q-grid.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre,expit

from frw_collisionless_mother import frozen_mother_background


class CollisionlessMotherFV:
    def __init__(self,nq=150,x_initial=.05,qmin=1e-5,qmax=40.,g1=4,gauss_n=16):
        if nq<5 or x_initial<=0 or qmin<=0 or qmax<=qmin or gauss_n<3:
            raise ValueError("Invalid comoving momentum discretization")
        self.nq=nq
        self.x_initial=x_initial
        self.g1=g1
        self.edges=np.r_[0.,np.geomspace(qmin,qmax,nq)]
        lo=self.edges[:-1]
        hi=self.edges[1:]
        nodes,weights=roots_legendre(gauss_n)
        self.q=lo[:,None]+.5*(hi-lo)[:,None]*(nodes[None,:]+1)
        self.weight=.5*(hi-lo)[:,None]*weights[None,:]
        self.q2_weight=self.q**2*self.weight
        self.volume=(hi**3-lo**3)/3
        # Average the genuine initial FD state against q²dq, rather
        # than sampling FD at shell mean and biasing even the number.
        initial=expit(-np.hypot(self.q,x_initial))
        self.f_i=np.sum(self.q2_weight*initial,axis=1)/self.volume
        self.number=(g1/(2*math.pi**2)) * np.dot(self.volume,self.f_i)
        if not np.all((self.f_i>=0)&(self.f_i<=1)):
            raise ArithmeticError("Initial quantum occupation invalid")

    def moments(self,x):
        if x<=0:
            raise ValueError("x=M1/T must be positive")
        e=np.hypot(self.q,x)
        scalar=self.q2_weight*self.f_i[:,None]
        pref=self.g1/(2*math.pi**2)
        energy=float(pref*np.sum(scalar*e))
        pressure=float(pref/3*np.sum(scalar*self.q**2/e))
        derivative=float(pref*np.sum(scalar*x/e))
        identity=x*derivative-energy+3*pressure
        physical_slope=-4+x*derivative/energy
        return {
           "x_M1_over_T":x,
           "n_over_T3":float(self.number),
           "rho_over_T4":energy,
           "pressure_over_T4":pressure,
           "w":pressure/energy,
           "analytic_d_rho_over_T4_dx":derivative,
           "x_d_rho_dx_minus_rho_plus_3P":identity,
           "d_ln_physical_rho_by_d_ln_a":physical_slope,
           "expected_FRW_energy_slope":-3*(1+pressure/energy),
           "comoving_occupation_is_constant":True,
        }

    def reference_error(self,x):
        numerical=self.moments(x)
        exact=frozen_mother_background(x,x_dec=self.x_initial,g1=self.g1)
        return {
            "n_relative_error":float(
                self.number/exact["comoving_mother_number_over_T3"]-1),
            "rho_relative_error":float(
                numerical["rho_over_T4"]/exact["dimensionless_rho_over_T4"]-1),
            "pressure_relative_error":float(
                numerical["pressure_over_T4"]/exact["dimensionless_pressure_over_T4"]-1),
        }


def run():
    cases=[]
    for nq in (30,50,100,200,400,600):
        model=CollisionlessMotherFV(nq=nq)
        snapshots=[model.moments(x) for x in (.05,1.,5.,20.)]
        for row in snapshots:
            assert abs(row["x_d_rho_dx_minus_rho_plus_3P"])<1e-10
            assert abs(row["d_ln_physical_rho_by_d_ln_a"]-
                       row["expected_FRW_energy_slope"])<1e-10
        cases.append({
            "comoving_q_bins":nq,
            "number_at_all_x":float(model.number),
            "continuum_relative_errors_at_x5":model.reference_error(5.),
            "snapshots":snapshots
        })
    return {
      "scope":"Exact finite-volume collisionless FRW number and continuity equations; independent continuum energy and pressure errors",
      "x_initial":.05,
      "q_bin_representation":"fixed q=p/T, phase-volume mean FD occupancy in each shell",
      "normalization":"g1/(2pi²) q²dq with g1=4",
      "status":"SOURCE-FREE MASSIVE MOTHER FRW BENCHMARK — NOT COUPLED BOSE 1<->2 FREEZEOUT",
      "cases":cases,
      "main_finding":"Bin-number normalization and physical FRW energy continuity can be exact while energy and pressure moments converge only with q-bin refinement.",
      "next":"Combine with physical on-shell 1<->2 collisions while preserving positive Pauli/Bose states, both comoving charges and the full FRW work/energy budget."
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/mother_fv_frw_moments.json")
    opts=p.parse_args()
    r=run()
    dest=Path(opts.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(r,indent=2),encoding="utf-8")
    print(json.dumps(r,indent=2))
