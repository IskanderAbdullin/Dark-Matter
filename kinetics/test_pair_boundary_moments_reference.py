"""Controlled tests of edge bookkeeping for synthetic signed pair events."""
import json
import math
import numpy as np
from massive_shell_moments import shell_weights
from pair_boundary_moments_reference import project_with_edge_moments


def verify():
    seed, n = 261126, 50000
    rng = np.random.default_rng(seed)
    m, b = 0.08, 0.18
    p = np.exp(rng.uniform(np.log(.01), np.log(25), n))
    c = rng.uniform(-1, 1, n)
    # q per daughter is always below 26 on this chosen event sample.
    q = np.r_[0, np.geomspace(1e-5, 26, 48)]
    p0 = np.sqrt((1-(m+b)**2)*(1-(m-b)**2))/2
    e = np.hypot(p, 1)
    xs = [e*(1+m*m-b*b)/2+p*p0*c, e*(1+b*b-m*m)/2-p*p0*c]
    out = {}
    for label, w in [("positive", np.ones(n)),
                     ("signed", np.where(rng.integers(0, 2, n), 1., -1.))]:
        numbers, energies, edge_counts = [], [], []
        for x, mass in zip(xs, (m,b)):
            nn, ee, _, _ = shell_weights(q, mass)
            means = ee/nn
            cnt, res = project_with_edge_moments(x, w, means, [mass, np.hypot(q[-1], mass)])
            numbers.append(math.fsum(cnt)+math.fsum(res[:,0]))
            energies.append(math.fsum(cnt*means)+math.fsum(res[:,1]))
            edge_counts.append(int(np.count_nonzero((x < means[0]) | (x > means[-1]))))
        norm_n, norm_e = math.fsum(abs(w)), math.fsum(abs(w*e))
        en_rel = abs(math.fsum(energies)-math.fsum(w*e))/norm_e
        nn_rel = max(abs(v-math.fsum(w))/norm_n for v in numbers)
        assert sum(edge_counts) == 30 and en_rel < 5e-15 and nn_rel < 5e-15
        out[label] = {"edge_events_by_species":edge_counts,
                      "pair_energy_residual_over_absolute_energy_flux":en_rel,
                      "max_number_residual_over_absolute_number_flux":nn_rel}
    # A low edge pair may have exactly zero signed net count but nonzero energy.
    counts, edges = project_with_edge_moments([1.01, 1.02], [1., -1.], [1.1, 1.4], [1., 2.])
    assert np.all(counts==0) and edges[0,0]==0 and abs(edges[0,1]+.01) < 2e-16
    try:
        project_with_edge_moments([2.01], [1.], [1.1, 1.4], [1., 2.])
    except ValueError as e:
        assert "physical grid" in str(e)
    else:
        raise AssertionError("physical overflow must not be silently accepted")
    return {"seed":seed, "events":n, **out, "signed_zero_count_nonzero_energy":True,
            "physical_overflow_rejected":True}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
