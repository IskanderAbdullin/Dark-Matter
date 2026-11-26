"""Validate the passive remap and match the archived instantaneous rate."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import points, snapshot_functions
from three_species_decay_onsager_reference import m, mp
from warm_ballistic_heavy_exposure_reference import redshift


def read(path):
    return json.loads(Path(path).read_bytes())


def check(exposure, fixed, paths):
    errors = []
    for path in paths:
        state = read(path)["snapshots"][-1]
        p = np.asarray(state["lab_momentum_edges_GeV"])
        f, _, g, _ = snapshot_functions(state)
        for b in (1.12, 1.57, 2.8, 3.11, 7.):
            for mass, fn in ((mp, f), (m, g)):
                p0, _, d0 = points(p, mass, fn, 24)
                p1, _, d1 = points(p/b, mass, redshift(fn, mass, b), 24)
                errors += [abs(b**3 * d1.sum()/d0.sum() - 1.),
                           abs(b**5 * (d1*p1*p1).sum()/(d0*p0*p0).sum() - 1.)]
    u = read(exposure)["cases"]
    v = read(fixed)["cases"]
    anchor_errors = []
    for a, c in zip(u, v, strict=True):
        old = c["orders"][-1]
        old_rates = np.asarray([old["forward_phi_phi_heavy_count_over_total_NH"],
                                old["forward_opposite_charge_chi_heavy_count_over_total_NH"]])
        old_rates *= a["retained_number_at_anchor_GeV3"] / a["normalization_at_anchor_GeV3"]
        actual = np.asarray(a["controls"][-1]["anchor_phi_and_opposite_stable_per_log_b"])
        np.testing.assert_allclose(actual, old_rates, rtol=.0008, atol=0.)
        anchor_errors.append(float(np.max(abs(actual/old_rates - 1.))))
    assert max(errors) < 2e-10
    return dict(max_passive_moment_relative_error=float(max(errors)),
                max_anchor_relative_errors=anchor_errors)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("exposure")
    p.add_argument("fixed")
    p.add_argument("paths", nargs=2)
    args = p.parse_args()
    print(json.dumps(check(args.exposure, args.fixed, args.paths), indent=2))
