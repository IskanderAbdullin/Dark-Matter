"""Independent coordinate and neighboring-shell controls for the cold reverse lift."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import expit, roots_legendre

from archived_scalar_pair_reaction_reference import gross, m, mp
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from warm_cold_light_finite_pair_reference import tagged
from warm_cold_light_pair_partition_reference import lifted
from warm_light_occupation_stage_audit_reference import digest, read


def synthetic():
    p = np.r_[0., .06, .11, .175, .27, .41, .65, .9]
    fe, be = np.hypot(p, m), np.hypot(p, mp)
    F = lambda e: np.where(e < fe[-1], expit((-.025 - e) / .017), 0.)
    B = lambda e: np.where(e < be[-1], 1. / np.expm1((e + .025) / .017), 0.)
    out = []
    for a, b in (((4,12,16,8,12),(5,16,10,6)), ((6,20,22,12,16),(8,26,12,8))):
        project, both, closure = lifted(p, F, B, fe, be, a)
        _, _, allp = gross(p, mp, m, B, F, True, *b, project=True)
        _, _, other = gross(p, mp, m, B, F, True, *b, project=True, min_final_energy=fe[1])
        subset = allp - other
        L1 = np.abs(project - subset).sum(axis=2) / project.sum(axis=2)
        frac = both / (project[0, 0].sum() / 2)
        assert np.max(L1) < .006 and .45 < frac < .55
        out.append(dict(orders=[a,b], relative_projection_L1=L1.tolist(),
                        native_tag_number_defect=float(subset[1,0,1]/project[1,0,1]-1.),
                        rate_both_cold_fraction=float(frac), lift_exact_tag_closure=closure))
    return out


def scan(path):
    row = read(path)
    row = row if isinstance(row, dict) else row[-1]
    S = row["snapshot"]
    assert row["snapshot_sha256"] == digest(S)
    p = np.asarray(S["lab_momentum_edges_GeV"])
    B, be, F, fe = snapshot_functions(S)
    g = []
    t, w = roots_legendre(4)
    for i in range(12):
        x = (p[i]+p[i+1])/2 + (p[i+1]-p[i])*t/2
        weights = 4/(2*np.pi**2)*(p[i+1]-p[i])*w*x*x/2
        v = np.array([tagged(u,F,B,fe,be,(8,6,8,8))[0] for u in x])
        g.append(weights @ v)
    a = np.asarray(g)
    out = []
    _, _, L = gross(p,m,mp,F,B,False,3,7,6,2,project=True)
    for order in ((3,7,6,2),(4,10,10,4)):
        _, G, R = gross(p,mp,m,B,F,True,*order,project=True)
        out.append(dict(order=order, gain_ratio=(R[1,0,1:13]/a[:,0]).tolist(),
                        loss_ratio_at_evolution_order=(L[0,0,1:13]/a[:,1]).tolist(),
                        tagged_shell_gain_fraction_of_native_total=(a[:,0]/(2*G[1])).tolist()))
    return dict(T_GeV=S["T_visible_GeV"], sha256=digest(S), edges_GeV=p[:13].tolist(), comparison=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("files", nargs="+", type=Path)
    a = p.parse_args()
    print(json.dumps(dict(double_cold_control=synthetic(), neighboring_shells=[scan(f) for f in a.files]),
                     indent=2, allow_nan=False))
