"""Replay diagnostics on the complete accepted two-packet stress checkpoints."""
import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from warm_cohort_two_packet_stress_reference import mom, run
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest, read


def audit(eos, work, tapes, vacuum, hard):
    checked = run(eos, work, tapes, vacuum, hard, False)
    source = read(work / "source.json")["snapshots"][0]
    stress = read(work / "q64" / "checkpoint.json")
    blank = read(vacuum)
    a = stress["snapshot"]["scale_factor"] / source["scale_factor"]
    assert abs(a / 1.001 - 1.) < 5e-14
    assert stress["plan"]["outer_log_scale_grid"] == blank["plan"]["outer_log_scale_grid"]
    assert stress["plan"]["eos_sha256"] == blank["plan"]["eos_sha256"] == sha(eos)
    assert stress["plan"]["local_prompt_jacobian"] and blank["plan"]["local_prompt_jacobian"]
    for row in (stress, blank):
        assert row["snapshot_sha256"] == digest(row["snapshot"])
    np.testing.assert_allclose(stress["snapshot"]["lab_momentum_edges_GeV"],
                               blank["snapshot"]["lab_momentum_edges_GeV"], rtol=2e-14, atol=0.)
    new = import_restart(stress["snapshot"]["positive_exterior_restart"]).array()
    assert new.shape[1] == 0 or np.all(new[5] >= 0.)
    # Old cohorts were moved into the explicit state. Only new receipts
    # are added back in this conservation test, never the inherited source.
    number_before = mom(source)[0].sum()
    number_after = (mom(stress["snapshot"])[0].sum() + new[5].sum()) * a**3
    assert abs(number_after / number_before - 1.) < 2e-13
    ledgers = [{k: row[k] for k in ("T_visible_GeV", "energy_ledger", "count_ledger",
                                   "accepted_microsteps", "rejected_microsteps")}
               for row in stress["ledger"]]
    assert max(max(abs(row[k]) for k in ("energy_ledger", "count_ledger")) for row in ledgers) < 3e-13
    initial = np.array(checked["initial_high_inverse_and_tagged_energy_rate"])[-1]
    evolved = np.array(checked["terminal_high_inverse_and_tagged_energy_rate"])[-1]
    control = np.array(checked["vacuum_terminal_high_inverse_and_tagged_energy_rate"])[-1]
    return dict(hard_GeV=hard, matched_visible_T_GeV=stress["snapshot"]["T_visible_GeV"],
                expansion_factor=a, inclusive_total_count_closure=float(number_after/number_before-1.),
                segment_ledgers=ledgers, terminal_over_initial=(evolved/initial).tolist(),
                terminal_over_vacuum=(evolved/control).tolist(), results=checked)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("eos", "work", "vacuum"):
        p.add_argument(name, type=Path)
    p.add_argument("--tapes", type=Path, default=Path("results"))
    p.add_argument("--hard", type=float, default=.58)
    a = p.parse_args()
    print(json.dumps(audit(a.eos,a.work,a.tapes,a.vacuum,a.hard),indent=2,allow_nan=False))
