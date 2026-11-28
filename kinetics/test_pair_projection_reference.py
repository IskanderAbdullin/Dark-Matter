"""Two-body on-shell Monte Carlo reference; not a collision solver."""
import json
import numpy as np
from massive_shell_moments import shell_weights
from pair_projection_reference import bracket_energy, deposit


def verify():
    seed, n = 261126, 50000
    rng = np.random.default_rng(seed)
    m, b = 0.08, 0.18    # Parent mass set to 1.
    parent_p = np.exp(rng.uniform(np.log(.01), np.log(25), n))
    cos = rng.uniform(-1, 1, n)
    parent_e = np.hypot(parent_p, 1)
    e0, r0 = (1+m*m-b*b)/2, (1+b*b-m*m)/2
    p0 = np.sqrt((1-(m+b)**2)*(1-(m-b)**2))/2
    phi, ferm = parent_e*e0+parent_p*p0*cos, parent_e*r0-parent_p*p0*cos
    assert np.max(abs(phi+ferm-parent_e)/parent_e) < 5e-16
    tight = np.r_[0, np.geomspace(1e-5, 26, 48)]
    dropped = 0
    for e, mass in ((phi, m), (ferm, b)):
        nn, ee, _, _ = shell_weights(tight, mass)
        dropped += int(np.count_nonzero((e < ee[0]/nn[0]) | (e > ee[-1]/nn[-1])))
    assert dropped > 0
    report = []
    for cells in (48, 96, 192):
        q = np.r_[0, np.geomspace(1e-5, 36, cells)]
        v = []
        for e, mass in ((phi, m), (ferm, b)):
            number, energy, _, _ = shell_weights(q, mass)
            means = energy/number
            # All sampled events must be covered by the positive bracket.
            i, j, w = bracket_energy(e, means)
            assert np.all((0 <= w) & (w <= 1))
            population = deposit(e, means)
            assert abs(population.sum()-n)/n < 5e-15
            assert abs(population @ means - e.sum())/e.sum() < 5e-15
            one = np.searchsorted(q, np.sqrt((e-mass)*(e+mass)), side="right")-1
            assert np.all((0 <= one) & (one < cells))
            v.append(((1-w)*means[i]+w*means[j], means[one]))
        linear = (v[0][0]+v[1][0]-parent_e)/parent_e
        naive = (v[0][1]+v[1][1]-parent_e)/parent_e
        assert np.max(abs(linear)) < 1e-15
        report.append({"cells": cells, "max_relative_pair_error": float(np.max(abs(linear))),
                       "single_cell_relative_rms": float(np.sqrt(np.mean(naive**2))),
                       "single_cell_max_relative_error": float(np.max(abs(naive)))})
    number, energy, _, _ = shell_weights([0., .01, .03], b)
    try:
        bracket_energy([b], energy/number)
    except ValueError as e:
        assert "outside" in str(e)
    else:
        raise AssertionError("rest energy must lie outside mean support")
    assert report[-1]["single_cell_relative_rms"] < report[0]["single_cell_relative_rms"]/3
    return {"seed": seed, "events": n, "daughter_masses_parent_units": [m,b],
            "meshes": report, "near_rest_support_failure_reproduced": True,
            "daughter_events_outside_tight_48_cell_mean_support": dropped}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
