"""Isolate pair and elastic effects on the same finite-parent trajectory."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit

from joint_prompt_two_moment_evolution_reference import Cells
from three_species_decay_onsager_reference import m,mp

from warm_clock_path_error_profile_reference import audited
from warm_finite_parent_expansion_species_budget_reference import measure
from warm_light_occupation_stage_audit_reference import digest,read


MODES=("expand-finite","expand-finite-pair","expand-finite-elastic","expand-finite-all")


def entropy(s):
    p=np.asarray(s["lab_momentum_edges_GeV"])
    out=[]
    for name,mass,g in (("scalar",mp,1),("fermion",m,4)):
        c=Cells(p,mass,g)
        t=np.asarray(s["entropy_variable_"+name])
        v=t[0,:,None]+t[1,:,None]*c.u
        if g==1:
            assert np.min(v)>0.
            f=np.exp(-v)/(-np.expm1(-v))
            ent=v*f-np.log(-np.expm1(-v))
        else:
            f=expit(-v)
            ent=v*f+np.logaddexp(0.,-v)
        out.append(np.sum(c.vol[:,None]*c.w*ent))
    return np.asarray(out)


def audit(start,folders):
    r=read(start)
    S=r["snapshot"]
    assert r["hash"]==digest(S)
    C=r["plan"]["normalized_cut_GeV"]
    N0,K0,E0,W0=measure(S,C)
    rho=E0.sum()
    S_entropy=entropy(S)
    common=None
    entries=[]
    values=[]
    for folder,mode in zip(folders,MODES,strict=True):
        t=audited(folder)
        rows=t["rows"]
        f=rows[-1]
        assert t["mode"]==mode and t["input_sha256"]==r["hash"] and t["cut_GeV"]==C
        assert abs(t["W_initial_over_rho"]/(sum(W0)/rho)-1.)<3e-14
        case_key=(t["eos_sha256"],f["c"],f["time_J"],f["qsha256"],f["wsha256"],[j["time_J"] for j in rows])
        if common is None:
            common=case_key
        assert common==case_key
        output=read(folder/f"step-{len(rows)}.json")["daughter_snapshot"]
        N,K,E,W=measure(output,C)
        a=f["c"]
        assert abs(W.sum()/rho/f["W_over_initial_daughter_rho"]-1.)<3e-14
        now=np.stack([a**3*N/N0,a**4*K/rho,a**4*E/rho,W/rho,a**3*entropy(output)/sum(S_entropy)])
        initial=np.stack([np.ones(2),K0/rho,E0/rho,W0/rho,S_entropy/sum(S_entropy)])
        assert ("pair_comoving_count_exchange_over_initial" in f)==(mode.endswith("-pair") or mode.endswith("-all"))
        assert ("elastic_panels" in f)==(mode.endswith("-elastic") or mode.endswith("-all"))
        defect=max(abs(f["relative_energy_defect"]),*map(abs,f["relative_N_defects"]))
        assert defect<3e-10
        entries.append(dict(mode=mode,phase_sha256=f["daughter_hash"],final_W_over_initial_rho=float(sum(W)/rho),
                            component_changes=(now-initial).tolist(),closure_defect=defect))
        values.append(now)
    X=np.array(values)
    contrasts=X[1:]-X[0]
    interaction=contrasts[2]-contrasts[0]-contrasts[1]
    # K already has the common normalization a^4/rho0; count is
    # per species normalized by its own initial comoving N. Keeping
    # the algebra in these quantities avoids an arbitrary energy unit
    # for the effective comoving kinetic energy of one particle.
    mean=X[:,1]/X[:,0]
    num_source=(X[1:,0]-X[0,0])*mean[0]
    per_particle_source=X[1:,0]*(mean[1:]-mean[0])
    np.testing.assert_allclose(num_source+per_particle_source,contrasts[:,1],
                               rtol=4e-7,atol=6e-16)
    return dict(scope="Matched short numerical block budget. The interaction is a discrete algebraic contrast, not a resolved mixed reaction rate.",
                hard_GeV=r["plan"]["hard_GeV"],input_hash=r["hash"],clock=common,
                columns=["scalar","light_fermion"],entries=entries,
                rows=["comoving_count_over_initial_per_species","comoving_kinetic_energy_over_initial_total",
                      "comoving_full_energy_over_initial_total","physical_W_over_initial_total","comoving_daughter_entropy_over_initial_daughter_entropy"],
                pair_only_minus_prompt=contrasts[0].tolist(),
                elastic_only_minus_prompt=contrasts[1].tolist(),
                both_minus_prompt=contrasts[2].tolist(),nonadditive_interaction=interaction.tolist(),
                kinetic_energy_contrast_decomposition=dict(
                   cases=["pair","elastic","both"],count_at_prompt_mean=num_source.tolist(),
                   changed_mean_at_new_count=per_particle_source.tolist(),
                   mean_per_particle_relative_change=(mean[1:]/mean[0]-1.).tolist()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("state",type=Path)
    p.add_argument("folders",type=Path,nargs=4)
    a=p.parse_args()
    print(json.dumps(audit(a.state,a.folders),indent=2,allow_nan=False))
