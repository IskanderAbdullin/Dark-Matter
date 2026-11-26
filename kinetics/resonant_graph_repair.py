"""Construct an invariant-clean adaptive resonant 1<->2 graph.

Problem: separated narrow daughter-energy bands create a fake modular
pseudoenergy, even with exact discrete two-number/energy conservation.

Construct a connected "fine infrared core + coarse energy comb": indices
0..(core-1) and multiples of stride from stride to max_index. The overlap
of sumsets is deliberately sufficient to constrain the collision-invariant
nullspace to exactly three *on the tested examples*; NEVER assume this for
other parameters. Check rank and parent-thermal moments independently.

This is a collision-network and parent-quadrature design tool, NOT a
kinetic freeze-out solver. Mother sums can exceed 10,000 states, so dense
BDF integration is intentionally NOT attempted or advertised.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from collision_invariant_rank import audit_integer_grid

DEFAULT_MASS = .1
DEFAULT_DELTA = .001
DEFAULT_CORE = 100
DEFAULT_STRIDE = 50


def bridged_indices(core=DEFAULT_CORE, stride=DEFAULT_STRIDE,
                    max_index=8500):
    if not all(isinstance(k, int) for k in (core, stride, max_index)):
        raise TypeError("all mesh indices must be integers")
    if core < 3 or stride < 2 or max_index <= stride:
        raise ValueError("invalid fine core / coarse stride / max index")
    return sorted(set(range(core)).union(range(stride, max_index+1, stride)))


def exact_mother_thermal_moments(mass=DEFAULT_MASS, g1=4):
    if mass <= 0 or g1 <= 0:
        raise ValueError("positive mass/degeneracy required")
    pref = g1/(2*math.pi**2)
    def thermal(p):
        E = math.hypot(p,mass)
        return expit(-E)
    n = pref*quad(lambda p: p*p*thermal(p),0,np.inf,epsabs=1e-12)[0]
    rho = pref*quad(lambda p: p*p*math.hypot(p,mass)*thermal(p),
                    0,np.inf,epsabs=1e-12)[0]
    pressure = (pref/3)*quad(
        lambda p:p**4/math.hypot(p,mass)*thermal(p),
        0,np.inf,epsabs=1e-12)[0]
    return {"number_over_T3":float(n),
            "energy_over_T4":float(rho),
            "pressure_over_T4":float(pressure)}


def parent_grid_moments(indices, delta=DEFAULT_DELTA,
                        mass=DEFAULT_MASS, g1=4):
    idx=np.array(sorted(set(int(v) for v in indices)),dtype=np.int64)
    if len(idx)<3 or np.any(idx<0) or delta<=0 or mass<=0 or g1<=0:
        raise ValueError("invalid grid or physical parameters")
    daughter=(idx+.5)*delta
    on_shell=4*daughter[:,None]*daughter[None,:]>=mass*mass
    mother_keys=np.unique((idx[:,None]+idx[None,:]+1)[on_shell])
    if len(mother_keys)<2:
        raise ValueError("too few allowed mother energies")
    E=mother_keys*delta
    first=max(mass,E[0]-(E[1]-E[0])/2)
    edges=np.r_[first,(E[:-1]+E[1:])/2,
                E[-1]+(E[-1]-E[-2])/2]
    p3=lambda x:np.maximum(x*x-mass*mass,0)**1.5
    volume=(g1/(6*math.pi**2))*(p3(edges[1:])-p3(edges[:-1]))
    if not np.all(volume>0):
        raise ValueError("nonpositive mother phase volume")
    f=expit(-E)
    number=float(np.dot(volume,f))
    energy=float(np.dot(volume*E,f))
    pressure=float(np.dot(volume*((E*E-mass*mass)/(3*E)),f))
    return {"active_mother_modes":int(len(E)),
            "parent_energy_min":float(E.min()),
            "parent_energy_max":float(E.max()),
            "max_mother_energy_step":float(np.max(np.diff(E))),
            "minimum_daughter_energy":float(daughter.min()),
            "maximum_daughter_energy":float(daughter.max()),
            "number_over_T3":number,"energy_over_T4":energy,
            "pressure_over_T4":pressure,
            "parent_shells_all_positive":True}


def evaluate_mesh(core=DEFAULT_CORE,stride=DEFAULT_STRIDE,max_index=8500,
                  delta=DEFAULT_DELTA,mass=DEFAULT_MASS,reference=None):
    idx=bridged_indices(core,stride,max_index)
    graph=audit_integer_grid(idx,delta=delta,mass=mass)
    if graph["analytical_three_invariant_residual"]>=1e-10:
        raise ArithmeticError("physical energy invariant itself failed")
    mom=parent_grid_moments(idx,delta=delta,mass=mass)
    ref=reference if reference is not None else exact_mother_thermal_moments(mass)
    errors={key:float(mom[key]/ref[key]-1)
            for key in ("number_over_T3","energy_over_T4","pressure_over_T4")}
    graph_is_clean=bool(graph["accept_for_unique_chemical_equilibrium"])
    return {
        "mesh":{"core":core,"stride":stride,"max_index":max_index,
                "delta":delta,"M1_over_T":mass,"daughter_modes":len(idx),
                "index_construction":"range(core) union range(stride,max_index+1,stride)"},
        "collision_invariant_nullity":graph["active_graph_collision_invariants"],
        "extra_collision_invariants":graph["additional_invariants_beyond_number_number_energy"],
        "active_reactions":graph["active_pair_reactions"],
        "isolated_daughters":graph["isolated_boson_modes"]+
                             graph["isolated_chi2_modes"],
        "graph_is_clean":graph_is_clean,
        "discrete_mother_moments":mom,
        "independent_continuum_reference":ref,
        "mother_relative_errors":errors,
        "parent_moment_max_abs_relative_error":max(abs(x) for x in errors.values()),
        "numerical_graph_and_moments_ready_for_further_studies":
            graph_is_clean and max(abs(x) for x in errors.values())<1e-3,
        "NOT_a_validated_cosmological_PSD":True,
        "remaining":"No demonstrated full nonlinear evolution, IR spectral convergence, zero mode, FRW, kinetic source Lagrangian or observational likelihood."
    }


def run(quick=False):
    qmax=(5000,8500) if quick else (5000,7000,8500,10000)
    ref=exact_mother_thermal_moments()
    results=[evaluate_mesh(max_index=m,reference=ref) for m in qmax]
    for r in results:
        assert r["graph_is_clean"],r
    return {
        "status":"INVARIANT-CLEAN CONNECTED ADAPTIVE GRAPH AND STATIC FD-MOMENT TEST; not a freeze-out spectrum",
        "baseline":"ea6ec461997600d5f00bf244ea808607d40ef8dd",
        "result":"The spurious-invariant obstruction is removable by bridging sumset connectivity, at the price of thousands of unique mother energies.",
        "cases":results,
        "important_caveat":"Number/moment accuracy at initial FD equilibrium is necessary but not sufficient for correct nonlinear collision rates. Solver cost and IR Bose occupation still require validation.",
    }


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    parser.add_argument("--output",default="outputs/resonant_graph_repair.json")
    opts=parser.parse_args()
    output=run(quick=opts.quick)
    path=Path(opts.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(output,indent=2),encoding="utf-8")
    print(json.dumps(output,indent=2))
