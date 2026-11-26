"""Audit a conservative reverse-coordinate cut across multiple cold shells."""

import argparse
import json
from pathlib import Path

import numpy as np

from frw_exterior_restart_reference import import_restart
from joint_prompt_two_moment_evolution_reference import Cells
from late_radiation_background_reference import MPL, TabulatedRadiation
from scalar_pair_expanding_continuation_reference import rate
from three_species_decay_onsager_reference import m
from warm_cold_light_pair_partition_reference import split_rate
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_reentry_shape_audit_reference import get_mom


def measure(eos, file, finer=False):
    entry = read(file)
    entry = entry[-1] if isinstance(entry, list) else entry
    S = entry["snapshot"]
    assert digest(S) == entry["snapshot_sha256"]
    p = np.array(S["lab_momentum_edges_GeV"])
    th = [S["entropy_variable_scalar"], S["entropy_variable_fermion"]]
    v = Cells(p,m,4).vol
    rho = TabulatedRadiation(eos, S["T_visible_GeV"]).state(1.)[1]
    grav = get_mom(p, th)[2].sum() + import_restart(S["positive_exterior_restart"]).energy(1.)
    H = np.sqrt((rho + grav)/3)/MPL
    native_o = (3,7,6,2)
    raw, _, _ = rate(p,th,native_o)
    normal = np.sum(abs(raw), axis=2)
    out = []
    methods = [(4,8,6,8,8), (6,12,9,12,12)]
    if finer:
        methods.append((9,18,14,18,16))
    for cut in (1,6,9):
        saved = []
        chosen = []
        for method in methods:
            answer,_,_,d = split_rate(p,th,native_o,method,cut_cells=cut)
            saved.append(answer)
            chosen.append(dict(tag_order=method, **d,
                normalized_rate_L1= (np.sum(abs(answer-raw),axis=2)/normal).tolist(),
                frozen_first_12_signed_per_H=(answer[1,0,1:13]/(H*v[:12])).tolist()))
        diffs = [np.sum(abs(v-u),axis=2)/normal for u,v in zip(saved[:-1],saved[1:],strict=True)]
        rates = [z[1,0,1:cut+1] for z in saved]
        max_relative = [max(abs(v-u)/np.maximum(np.abs(v),1.e-280))
                        for u,v in zip(rates[:-1],rates[1:],strict=True)]
        out.append(dict(cut_cells=cut, cut_GeV=p[cut], quadrature_change_L1=[x.tolist() for x in diffs],
                        worst_tagged_cell_signed_relative_quadrature_change=[float(x) for x in max_relative],
                        controls=chosen))
    return dict(name=file.name, T_GeV=S["T_visible_GeV"], sha256=digest(S), p_edges_GeV=p[:13].tolist(),
                native_first_12_signed_per_H=(raw[1,0,1:13]/(H*v[:12])).tolist(), cuts=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos", type=Path)
    p.add_argument("files", nargs="+", type=Path)
    p.add_argument("--finer", action="store_true")
    args = p.parse_args()
    print(json.dumps([measure(args.eos, f, args.finer) for f in args.files], indent=2, allow_nan=False))
