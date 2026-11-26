"""Admit all recorded exterior cohorts into the isolated warm recycling bound."""
import argparse
from pathlib import Path
import json

import numpy as np

from frw_exterior_restart_reference import import_restart
from three_species_decay_onsager_reference import M
from warm_archived_expansion_exterior_envelope_reference import MASSES, from_ledger, history
from warm_early_support_time_interaction_reference import read, sha
from warm_early_serial_support_audit_reference import rows
from warm_isolated_recycling_energy_envelope_reference import moments, radiation_bound


def run(eos, cores, checkpoints, seeds, producers, budget, isolated, old_audit):
    source, closed, empirical = [read(p) for p in (budget, isolated, old_audit)]
    goal = closed["target_total"]
    K = closed["combined_K_GeVminus2"]
    ratio = MASSES[0]/MASSES[1]
    assert empirical["EOS_sha256"] == sha(eos)
    out = []
    for i, (core, cp, seed, prod) in enumerate(zip(cores, checkpoints, seeds, producers, strict=True)):
        q = (16, 32)[i]
        C = read(cp)
        S = C["snapshot"]
        init = read(seed)["snapshots"][-1]
        T = S["T_visible_GeV"]
        a0, a = init["scale_factor"], S["scale_factor"]
        assert S == read(core)["snapshots"][-1]
        assert C["blocks_total"] == C["blocks_finished"] == 16
        assert C["plan"]["eos_sha256"] == sha(eos)
        p = read(prod)["rows"][-1]
        assert abs(p["final_visible_temperature_GeV"]/init["T_visible_GeV"]-1.) < 1e-12
        initial_n = np.asarray(init["represented_species_count_a3"])
        norm = np.asarray(source["cases"][i]["retained_fixed_Omega_h2"])/(MASSES*initial_n)
        old_omega = np.asarray(source["cases"][i]["exterior_fixed_Omega_h2"])
        assert abs(norm[0]/norm[1]-1.) < 3e-12
        old_number = old_omega/(MASSES*norm)
        old = from_ledger(p["visible_bath_feedback"], a0**4, old_number, a0)
        old_ballistic, _ = history([old], a)
        # If interactions redistribute old momentum only within a separate
        # pre-anchor reservoir, positivity of pressure is enough instead.
        old_pressure_only_max = (old[1]+old[3])*(a0/a)**3
        assert old_ballistic[1] <= old_pressure_only_max
        V = empirical["cases"][i]
        assert V["initial_support"] == q
        assert abs(V["data"][-1]["photon_T_GeV"]/T-1.) < 1e-12
        np.testing.assert_allclose(old_ballistic, V["data"][-1]["old_ext_rho_interval_GeV4"],
                                   rtol=2e-12, atol=0.)
        W = rows(init, C, source["cases"][i], sha(seed))[-1]
        Y = (np.array([1., ratio]) @ W["inclusive_fixed_count_omega"])
        prior = closed["cases"][i]
        np.testing.assert_allclose(Y, prior["conservative_initial_Y"], rtol=2e-10, atol=0.)
        assert prior["q"] == q
        N, rho = moments(S).sum(axis=0)
        tails = import_restart(S["positive_exterior_restart"])
        new_count_fraction = float(np.sum(tails.array()[5])/N)
        new_upper = tails.upper_bound(1.)
        clock = radiation_bound(eos, T, 500.)
        time = clock["upper_integral_times_Hph"]/clock["photon_anchor_GeV"]
        gain = []
        for label, extra in (("retained_only", 0.),
                             ("all_known_old_ballistic_and_new", old_ballistic[1]+new_upper),
                             ("all_known_old_pressure_only_and_new", old_pressure_only_max+new_upper)):
            # N is deliberately just the core normalizer. Multiply by the
            # larger inclusive Y to cap absolute gains conservatively.
            multiplier = 1 + K*(rho+extra)**2/(M*M*N)*time
            Omega = (Y*multiplier-ratio*goal)/(1-ratio)
            gain.append(dict(history=label, added_energy_over_core=float(extra/rho),
                             inventory_ratio_upper=float(multiplier),
                             scalar_Omega_h2_upper_if_total_fixed=float(Omega)))
        control = prior["bounds"][1]["upper_inventory_multiplier"]
        assert abs(gain[0]["inventory_ratio_upper"]/control-1.) < 2e-14
        out.append(dict(q=q, anchor_T_GeV=T, known_new_count_over_core=new_count_fraction,
                         known_old_count_over_core=float(sum(old_number)/(a**3*N)),
                         new_exterior_rho_upper_GeV4=new_upper,
                         old_exterior_rho_ballistic_interval_GeV4=old_ballistic.tolist(),
                         old_exterior_rho_pressure_only_upper_GeV4=old_pressure_only_max,
                         same_total=goal, bounds=gain))
    return dict(scope="After the anchor, known old and new cohorts can re-enter "
                      "and interact arbitrarily inside the specified isolated network. "
                      "Does not redo evolution or allow transfer from the core into old "
                      "tails before the anchor, new visible energy, or unknown sources.",
                input_sha256=[sha(p) for p in [eos,*cores,*checkpoints,*seeds,*producers,budget,isolated,old_audit]],
                cases=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("eos", "budget", "isolated", "old_audit"):
        p.add_argument(key, type=Path)
    for key in ("cores", "checkpoints", "seeds", "producers"):
        p.add_argument("--"+key, nargs=2, required=True, type=Path)
    a = p.parse_args()
    print(json.dumps(run(a.eos, a.cores, a.checkpoints, a.seeds, a.producers,
                         a.budget, a.isolated, a.old_audit), indent=2, allow_nan=False))
