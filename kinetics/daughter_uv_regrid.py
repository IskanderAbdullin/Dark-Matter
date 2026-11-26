"""Number-exact positive Bose/Pauli regridding across expanding daughter UV grids.

A conserved massless daughter momentum is q=p/T when gstarS=const.
When the collision graph requires a larger daughter qmax at larger x,
a finite-volume regrid can preserve EVERY species' q²dq particle count
while inserting VACUUM in the newly represented high-q domain.

Bose f>=0 has NO Pauli cap; Fermi 0<=f<=1. Both are convexly remapped
with shell overlap in u=q³; counts are exact, but energy/entropy and
the momentum spectrum are altered by coarse graining. Critically,
zero occupation in added UV shells is an ASSUMPTION, not a recovered
physical thermal tail. This is a structural proof of topology-change
number conservation, not a cosmological rethermalization prescription.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import expit

from parent_q_remap import moving_parent_edges, energy
from sparse_kinetic_jacobian import SparseClosedDecay
from resonant_graph_repair import bridged_indices


def positive_species_regrid(old_edges,new_edges,occupation,statistics="bose",g=1):
    a=np.asarray(old_edges,dtype=float)
    b=np.asarray(new_edges,dtype=float)
    f=np.asarray(occupation,dtype=float)
    if statistics not in ("bose","fermi") or g<=0:
        raise ValueError("Bose or Fermi statistics and positive degeneracy required")
    if (len(a)<2 or len(b)<2 or a[0]!=0 or b[0]!=0
            or np.any(np.diff(a)<=0) or np.any(np.diff(b)<=0)
            or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b))):
        raise ValueError("Strictly increasing finite momentum grids starting at zero required")
    if b[-1]<a[-1]-1e-12:
        raise ValueError("Shrinking UV domain could delete physical particle number")
    if f.shape!=(len(a)-1,) or not np.all(np.isfinite(f)) or np.any(f<0):
        raise ValueError("Invalid positive species occupations")
    if statistics=="fermi" and np.any(f>1):
        raise ValueError("Pauli violation in the old distribution")
    ua=a**3;ub=b**3
    deposition=np.zeros(len(b)-1)
    i=j=0
    while i<len(f) and j<len(deposition):
        left=max(ua[i],ub[j])
        right=min(ua[i+1],ub[j+1])
        if right>left:deposition[j]+=f[i]*(right-left)
        if ua[i+1]<=ub[j+1]:i+=1
        else:j+=1
    newf=deposition/np.diff(ub)
    prior=float(np.dot(np.diff(ua),f))
    after=float(np.dot(np.diff(ub),newf))
    if abs(prior-after)>1e-11*max(1,prior):
        raise ArithmeticError("Species q³ particle number is not conserved")
    if min(newf)<-1e-12 or (statistics=="fermi" and max(newf)>1+1e-12):
        raise ArithmeticError("Quantum bounds violated")
    return newf


def number(qedges,occ,g=1):
    q=np.asarray(qedges,dtype=float)
    return float(g/(6*np.pi**2)*np.dot(np.diff(q**3),occ))


def massless_rho(qedges,occ,g=1):
    q=np.asarray(qedges,dtype=float)
    return float(g/(8*np.pi**2)*np.dot(np.diff(q**4),occ))


def parent_qedges(model,x):
    E=model.parent_E
    if x>=E[0]:
        raise ValueError("Energy center below physical rest mass")
    top=E[-1]+.5*(E[-1]-E[-2])
    return moving_parent_edges(E,x,math.sqrt(top*top-x*x)+1e-12)


def remap_three_species(old_model,new_model,old_state,xold,xnew):
    n=old_model.nq
    phi,chi,mother=old_model.unpack(np.asarray(old_state,dtype=float))
    if phi.shape!=(n,) or chi.shape!=(n,) or mother.shape!=(old_model.nparent,):
        raise ValueError("Three-species state does not match old on-shell model")
    if (old_model.graph["active_graph_collision_invariants"]!=3 or
            new_model.graph["active_graph_collision_invariants"]!=3):
        raise ValueError("Reaction graph has spurious invariants")
    oldq=parent_qedges(old_model,xold)
    newq=parent_qedges(new_model,xnew)
    if newq[-1]<oldq[-1]:
        raise ValueError("Parent momentum UV support shrank")
    # The collision ODE must use the SAME parent q³ phase volumes as
    # this Liouville remap; otherwise its conserved charges would be
    # normalized on a different set of truncated mother shells.
    old_model.vol_parent=4/(6*np.pi**2)*np.diff(oldq**3)
    new_model.vol_parent=4/(6*np.pi**2)*np.diff(newq**3)
    phi_new=positive_species_regrid(old_model.edges,new_model.edges,phi,"bose",g=1)
    chi_new=positive_species_regrid(old_model.edges,new_model.edges,chi,"fermi",g=4)
    mom_new=positive_species_regrid(oldq,newq,mother,"fermi",g=4)
    species_before={
        "phi":number(old_model.edges,phi,1),
        "chi2":number(old_model.edges,chi,4),
        "chi1":number(oldq,mother,4)}
    species_after={
        "phi":number(new_model.edges,phi_new,1),
        "chi2":number(new_model.edges,chi_new,4),
        "chi1":number(newq,mom_new,4)}
    charges_before=(species_before["chi1"]+species_before["phi"],
                    species_before["chi1"]+species_before["chi2"])
    charges_after=(species_after["chi1"]+species_after["phi"],
                   species_after["chi1"]+species_after["chi2"])
    charge_error=float(max(abs(a-b) for a,b in zip(charges_before,charges_after)))
    old_ode_charges=old_model.moments(old_state)
    new_state=np.r_[phi_new,chi_new,mom_new]
    new_ode_charges=new_model.moments(new_state)
    ode_charge_disagreement=max(
        abs(old_ode_charges["charge_chi1_plus_phi"]-charges_before[0]),
        abs(old_ode_charges["charge_chi1_plus_chi2"]-charges_before[1]),
        abs(new_ode_charges["charge_chi1_plus_phi"]-charges_after[0]),
        abs(new_ode_charges["charge_chi1_plus_chi2"]-charges_after[1]))
    if ode_charge_disagreement>1e-10:
        raise ArithmeticError("Collision ODE phase volumes disagree with Liouville charges")
    if charge_error>1e-10:
        raise ArithmeticError("Full-state remap lost closed reaction charge")
    def phys(old,b,c,p,q_edges,x):
        return massless_rho(old.edges,b,1)+massless_rho(old.edges,c,4)+energy(q_edges,p,x)
    old_energy_at_xnew=phys(old_model,phi,chi,mother,oldq,xnew)
    new_energy=phys(new_model,phi_new,chi_new,mom_new,newq,xnew)
    work=old_energy_at_xnew-phys(old_model,phi,chi,mother,oldq,xold)
    return {
       "remapped_state":new_state,
       "ode_moment_charge_consistency_residual":float(ode_charge_disagreement),
       "number_before":species_before,
       "number_after":species_after,
       "before_two_charges":charges_before,
       "after_two_charges":charges_after,
       "maximum_charge_drift":charge_error,
       "daughter_modes_before":old_model.nq,
       "daughter_modes_after":new_model.nq,
       "mother_modes_before":old_model.nparent,
       "mother_modes_after":new_model.nparent,
       "physical_free_expansion_energy_change":work,
       "pure_regridding_physical_energy_error":new_energy-old_energy_at_xnew,
       "added_momentum_domain_assumed_vacuum":True,
       "both_models_parent_phase_volumes_rebased_on_q3":True,
       "not_an_energy_conserving_remap":True,
       "not_CLASS_ready":True,
    }


def run():
    before=SparseClosedDecay(bridged_indices(max_index=300),M1_T=.02)
    after=SparseClosedDecay(bridged_indices(max_index=500),M1_T=.03)
    state=np.r_[2*np.exp(-before.e/.04),
                .1*np.exp(-before.e),
                expit(-before.parent_E)]
    result=remap_three_species(before,after,state,.02,.03)
    spectrum=result.pop("remapped_state")
    # After the topology change, test one genuinely CLOSED collision
    # substep with the NEW on-shell graph, same parent shell volumes.
    from scipy.integrate import solve_ivp
    previous=after.moments(spectrum)
    collision=solve_ivp(after.rhs_sparse,(0.,.0001),spectrum,method="BDF",
                        jac=after.jacobian,rtol=1e-7,atol=1e-12)
    if not collision.success:
        raise RuntimeError(collision.message)
    newstate=collision.y[:,-1]
    newmom=after.moments(newstate)
    keys=("charge_chi1_plus_phi","charge_chi1_plus_chi2","total_discrete_energy")
    errors={key:float(newmom[key]-previous[key]) for key in keys}
    if max(abs(v) for v in errors.values())>1e-8:
        raise ArithmeticError("New-grid collision after remap broke physical charges/discrete energy")
    fermi_chi=newstate[after.nq:2*after.nq]
    fermi_parent=newstate[2*after.nq:]
    if min(newstate)<-1e-8 or max(fermi_chi)>1+1e-8 or max(fermi_parent)>1+1e-8:
        raise ArithmeticError("New-grid collision violated quantum positivity")
    result["new_graph_post_remap_collision_test"]={
        "tGamma":.0001,
        "nfev":int(collision.nfev),
        "new_total_phi_number":newmom["n_phi"],
        "discrete_charge_and_energy_residuals":errors,
        "quantum_bounds_pass":True,
        "warning":"Short collision continuation only; newly introduced UV cells were vacuum and physical energy/entropy is not proven closed."
    }
    result["remapped_occupations_length"]=len(spectrum)
    result["max_boson_occupation"]=float(np.max(spectrum[:after.nq]))
    result["warning"]="Vacuum-added high-q support is a numerical modeling assumption; no nonlinear collision reinitialization, FRW energy closure or physical primordial spectrum."
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",default="outputs/three_species_regrid.json")
    args=p.parse_args()
    report=run()
    path=Path(args.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))
