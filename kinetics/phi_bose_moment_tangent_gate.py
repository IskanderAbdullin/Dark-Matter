"""Near-critical Bose inverse-fugacity moment feasible tangent-cone audit.

For physical Galerkin moments M_i =∫dΠ h_i b(Σ h_j η_j), J=dM/dη
is a symmetric NEGATIVE Gram matrix. If an IR hat has no collisional
source but its neighbor does, its inverse-J nodal eta derivative may
be negative: ηphi(rest) can hit zero even with ONLY net φ destruction.
This is a mathematical discretization/tangent obstruction, NOT BEC.
"""
from __future__ import annotations
import math
import numpy as np
from kinetics.phi_bose_moment_galerkin import Mode,physical_moments
from kinetics.phi_bose_nonuniform_ir_grid import nodes
from kinetics.phi_bose_implicit_moment_step import (
    physically_weighted_yukawa_events,source_with_jacobian,
)


def two_hat_example(a=3.,b=1.,c=2.,q=1.):
    """Exact Gram G>0, J=-G, source C=(0,-q)."""
    D=a*c-b*b
    if not(a>0 and b>0 and c>0 and D>0 and q>0):
        raise ValueError('SPD overlapping Bose hat Gram and positive removal q')
    G=np.array([[a,b],[b,c]],float)
    C=np.array([0.,-q])
    eta_rate=np.linalg.solve(-G,C)
    predicted=np.array([-b*q/D,a*q/D])
    return {'G':G.tolist(),'C':C.tolist(),
            'eta_rate':eta_rate.tolist(),'analytic_rate':predicted.tolist(),
            'rest_eta_decreases_despite_zero_rest_hat_source':bool(eta_rate[0]<0)}


def tangent(bins=16,scheme='power2',*,mphi=86e-9,pmax=2.,nr=12,nangle=12):
    upper=math.hypot(.1,pmax)+.2
    if mphi<=0 or mphi>=.08:raise ValueError('positive open Yukawa mass')
    modes=[Mode(m,nodes(m,upper,bins),degeneracy=4,fermion=True)
           for m in (.1,.02)]+[Mode(mphi,nodes(mphi,upper,bins,scheme))]
    eta=[modes[0].nodes/.9,modes[1].nodes,modes[2].nodes]
    events,_=physically_weighted_yukawa_events(mphi=mphi,pmax=pmax,nr=nr,nangle=nangle)
    C,_,H=source_with_jacobian(modes,eta,events)
    n=bins+1
    Cp=C[2*n:]
    moment=physical_moments(modes[2],eta[2],order=24)
    J=moment['dM_deta']
    eta_rate=np.linalg.solve(J,Cp)
    negative=[(i,float(eta[2][i]/(-eta_rate[i])))
              for i in range(n) if eta_rate[i]<0]
    idx,tlimit=min(negative,key=lambda z:z[1]) if negative else (-1,math.inf)
    energies=np.concatenate([mo.nodes for mo in modes])
    return {'bins':bins,'method':scheme,'mphi_GeV':mphi,
       'phi_rest_inverse_fugacity':eta[2][0],
       'phi_rest_collisional_moment_source':float(Cp[0]),
       'phi_total_collision_rate':float(np.sum(Cp)),
       'positive_scalar_hat_sources':int(np.count_nonzero(Cp>0)),
       'physical_event_entropy':float(H),
       'phi_rest_eta_time_derivative':float(eta_rate[0]),
       'first_three_phi_eta_derivatives':list(map(float,eta_rate[:3])),
       'linearized_boundary_limited_node':idx,
       'linearized_max_allowed_time':tlimit,
       'susceptibility_condition_number':moment['cond_number'],
       'total_energy_defect':float(np.dot(energies,C)),
       'dark_fermion_number_defect':float(sum(C[:n])+sum(C[n:2*n])),
       'parent_plus_phi_number_defect':float(sum(C[:n])+sum(C[2*n:])),
       'all_yukawa_net_phi_is_destruction_no_real_BEC':bool(np.all(Cp<=0)),
       'tangent_test_not_physical_cosmology_or_condensate':True}


def audit():
    rows=[tangent(n,scheme) for n in (8,16,32,64)
          for scheme in ('uniform','power2','power3','two_zone')]
    mass_scaling=[tangent(16,'power2',mphi=m) for m in (86e-9,86e-8,86e-7,86e-6)]
    return {'schema_version':1,'two_hat_exact_counterexample':two_hat_example(),
        'finite_tree_yukawa_physical_Gamma_but_synthetic_thermal_states':True,
        'rows':rows,'scalar_mass_scaling':mass_scaling,
        'no_convergence_BEC_or_FRW_claim':True}


if __name__=='__main__':
 import json
 print(json.dumps(audit(),indent=2))
