"""Audit endpoint counts across contiguous warm runs without duplicating tails."""
import argparse
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def endpoint(row, initial_retained, carry):
    n = np.asarray(initial_retained, float)
    old = np.asarray(carry, float)
    v = row["endpoint"]
    retained, inclusive = (n * np.asarray(v[k]) for k in
                            ("retained_counts_relative_to_start",
                             "inclusive_counts_relative_to_start"))
    new_tail = inclusive - retained
    assert np.all(retained > 0) and np.all(new_tail >= -2e-15 * n)
    defect = abs(inclusive.sum() / n.sum() - 1.)
    assert max(defect, abs(row["energy_ledger"]),
               abs(row["collision_number_ledger"])) < 5e-9
    return retained, old + new_tail


def audit(budget, start_shape, results, bridges):
    b = read(budget)[-1]
    start = read(start_shape)["snapshots"][-1]
    assert b["enabled"] and b["signed_visible_portal"]
    mass = np.array([86e-9, .02])
    n0 = np.asarray(start["represented_species_count_a3"])
    norm = np.asarray(b["retained_fixed_Omega_h2"]) / (mass * n0)
    assert abs(norm[0] / norm[1] - 1.) < 3e-12
    ex0 = np.asarray(b["exterior_fixed_Omega_h2"]) / (mass * norm)
    out = []
    for fname, shape in zip(results, bridges):
        r = read(fname)["rows"][-1]
        final = read(shape)["snapshots"][-1]
        n, ex = endpoint(r, n0, ex0)
        identities(r, final, start, n)
        out.append(row(r, n, ex, norm, mass))
    return dict(endpoint_cases=out, conditional_on_unchanged_visible_entropy=True)


def identities(r, final, start, n):
    assert np.isclose(r["initial_temperature_GeV"], start["T_visible_GeV"], rtol=1e-12, atol=0.)
    assert np.isclose(final["scale_factor"] / start["scale_factor"], r["endpoint"]["b"], rtol=1e-12, atol=0.)
    assert np.isclose(final["T_visible_GeV"], r["endpoint"]["photon_temperature_GeV"], rtol=1e-12, atol=0.)
    assert np.max(abs(np.asarray(final["represented_species_count_a3"])/n-1.)) < 3e-12


def row(r, n, ex, norm, mass):
    d = r["endpoint"]
    return dict(n_steps=r["steps"], temperature_GeV=d["photon_temperature_GeV"],
                retained_and_all_particle_tail_Omega_h2=((n + ex) * norm * mass).tolist(),
                all_tail_Omega_h2=(ex * norm * mass).tolist(),
                max_resolved_scalar_occupation=d["maximum_resolved_scalar_occupation"],
                scalar_fraction_of_own_saturation=d["scalar_count_over_own_saturation"],
                total_particle_count_identity_max_defect=max(abs(r["count_ledger"]),
                                                             abs(r["collision_number_ledger"])))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("budget")
    p.add_argument("start_shape")
    p.add_argument("--results", nargs="+", required=True)
    p.add_argument("--bridges", nargs="+", required=True)
    a = p.parse_args()
    assert len(a.results) == len(a.bridges)
    print(json.dumps(audit(a.budget, a.start_shape, a.results, a.bridges),
                     indent=2, allow_nan=False))
