"""Compare the native and lifted trajectories from identical warm restarts."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from warm_coupled_serial_restart_reference import sha
from warm_history_lower_matched_audit_reference import measures
from warm_light_occupation_stage_audit_reference import digest, read
from warm_near_rest_light_inverse_decay_reference import window
from three_species_decay_onsager_reference import mp


def tape(eos, seed, path, steps):
    S = read(seed)[-1]["snapshot"]
    prior = digest(S)
    a0 = S["scale_factor"]
    N0 = measures(S)[0].sum()
    old = import_restart(S["positive_exterior_restart"]).array()
    p = np.array(S["lab_momentum_edges_GeV"])
    hist = [measures(S)]
    obs = []
    start = []
    for i in range(steps+1):
        if i:
            row = read(path / f"step-{i}.json")
            assert row["index"] == i and row["previous"] == prior
            assert row["plan"]["tape_sha256"] == sha(seed)
            assert row["plan"]["EOS_sha256"] == sha(eos)
            S = row["snapshot"]
            prior = digest(S)
            assert prior == row["snapshot_sha256"]
        a = S["scale_factor"]/a0
        assert abs(a/1.10**(i/24)-1.) < 1e-13
        grid = np.array(S["lab_momentum_edges_GeV"])
        np.testing.assert_allclose(grid*a,p,rtol=1e-13,atol=0.)
        assert np.hypot(mp,grid[-1]) < window(grid[1])[-2]
        info = import_restart(S["positive_exterior_restart"]).array()
        n0 = old.shape[1]
        orig = info[:,:n0].copy()
        orig[[0,5]] *= np.array([a,a**3])[:,None]
        np.testing.assert_allclose(orig,old,rtol=1.5e-13,atol=0.)
        V = measures(S)
        closure = (V[0].sum()+info[5,n0:].sum()) * a**3/N0-1
        assert abs(closure)<5e-13
        occ = V[0,1,:12]/(4*np.diff(grid[:13]**3)/(6*np.pi**2))
        rec = dict(step=i, snapshot_sha256=prior, T_GeV=S["T_visible_GeV"], expansion=a,
                   inclusive_count_residual=float(closure), first_12_mean=occ.tolist())
        if i:
            rec.update(local_energy_defect=row["operator_summary"]["energy_ledger"],
                       slope=((occ-start)/np.log(a/obs[-1]["expansion"])).tolist(),
                       own_history=row["audit"]["solver"]["used_history"])
            hist.append(V*a**3)
        obs.append(rec)
        start = occ
    return obs,np.array(hist)


def compare(eos, seed, native, changed, steps):
    key = read(changed / "cold_partition_manifest.json")
    assert key["seed_sha256"] == sha(seed) and key["eos_sha256"] == sha(eos)
    a, X = tape(eos,seed,native,steps)
    b, Y = tape(eos,seed,changed,steps)
    np.testing.assert_allclose([x["T_GeV"] for x in a],[x["T_GeV"] for x in b],rtol=1e-13)
    N = X[0].sum(axis=2)
    L1 = abs(X-Y).sum(axis=3)/N
    summary = {}
    for tag, value in (("native",a),("lifted",b)):
        first = np.array([j["slope"][0] for j in value[1:]])
        summary[tag] = dict(first_cell_positive_steps=[x+1 for x in np.flatnonzero(first>0).tolist()],
                            cold_slope_range=[float(np.min(first)),float(np.max(first))],
                            max_count_defect=float(max(abs(j["inclusive_count_residual"]) for j in value)),
                            max_atomic_energy_defect=float(max(abs(j["local_energy_defect"]) for j in value[1:])))
    assert steps > 0
    return dict(scope="Matched selected warm comparison from common seeds; no physical abundance inference.",
                matched_steps=steps, manifest=key, native=a, lifted=b, summary=summary,
                L1_number_energy_pressure=L1.tolist(),
                max_L1=float(np.max(L1)))


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","seed","native","changed"):
        p.add_argument(key,type=Path)
    p.add_argument("--steps",type=int,required=True)
    args=p.parse_args()
    print(json.dumps(compare(args.eos,args.seed,args.native,args.changed,args.steps),
                     indent=2,allow_nan=False))
