"""Verify a cooling-history contrast at identical collisionless metric probes."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from export_finite_cooling_class_psd_reference import smooth
from export_finite_cooling_pair_mixture_reference import support
from finite_cooling_cdm_los_survey_reference import K, sha
from finite_cooling_reheat_portal_reference import TF
from three_species_decay_onsager_reference import mp

KEY = "los_over_cdm_joined_prompt_recoil_losIBP_direct_monopoleIC"


def recoil(eos, pairs, directory, sigma=.14):
    if sigma not in (.10,.14):
        raise ValueError("Match one of the published widths: 0.10 or 0.14.")
    record = json.loads(pairs.read_text())
    assert sha(eos) == record["eos_sha256"]
    bg = Background(eos)
    qfac = float(bg.scale(TF)) * TF / bg.T0
    factor = 2 * np.pi**2 / 45 * bg.gs0 / (mp * conversion(eos)["Omega_h2_per_GeV_times_yield"])
    directory.mkdir(parents=True, exist_ok=True)
    result = {}
    for ep in record["episodes"]:
        ratio = ep["Tmax_GeV"] / ep["T_R_GeV"]
        case = "instant" if ratio == 1. else "ratio5"
        assert abs(ratio - (1. if case == "instant" else 5.)) < 1e-12
        assert case not in result
        v = ep["modes"][-1]
        M = v["new_scalar_recoil_shape"]
        omega = float(np.sum(np.asarray(v["new_scalar_omega"])[:, 0]))
        p, w = support(np.sum(M["count_fraction_by_epoch"], axis=0),
                       np.sum(M["q2_weighted_count_fraction_by_epoch"], axis=0), factor * omega)
        q, f = smooth(qfac*p, w, sigma, 5200)
        f *= w.sum() / (4*np.pi*np.trapezoid(q*q*f, q))
        defect = float(4*np.pi*np.trapezoid(q**4*f,q)/(w @ (qfac*p)**2) - 1)
        assert np.all(f > 0.) and abs(defect) < 1e-7
        fs = directory / f"{case}-scalar-recoil-sigma{sigma:.2f}.dat"
        np.savetxt(fs, np.c_[q, f], fmt="%.16e")
        result[case] = dict(rest_omega_h2=omega, table_sha256=sha(fs), relative_q2_moment_defect=defect)
    assert len(result) == 2
    out = dict(components=result, source_sha256=sha(pairs), eos_sha256=sha(eos), q_factor=qfac)
    if sigma != .14:
        out["sigma_ln_q"] = sigma
    return out


def compare(finite_paths, instant_paths, finite_high, reconstruction, exports):
    rr, exp = (json.loads(p.read_text()) for p in (reconstruction, exports))
    source = {}
    for case in ("ratio5", "instant"):
        c = next(r for r in exp["components"] if r["case"] == case and r["sigma_ln_q"] == .14)
        h = rr["components"][case]
        assert abs(c["first_order_pair_omega_h2"]/h["rest_omega_h2"] - 1) < 1e-13
        source[case] = dict(rest_omega_h2=c["actual_fixed_count_Omega_h2"], recoil_omega_h2=h["rest_omega_h2"],
                            anchor_p_over_T_rms=c["reconstruction_anchor_p_over_T_rms"],
                            joined_sha256=c["sha256"], recoil_sha256=h["table_sha256"])
    def get(path, name):
        d = json.loads(path.read_text())
        r = d["results"]
        assert d["eos_sha256"] == rr["eos_sha256"] == exp["eos_sha256"]
        assert d["spectrum_sha256"][::2] == [source[name]["joined_sha256"], source[name]["recoil_sha256"]]
        C = np.array([[u[0] for u in x["cases"][-1][KEY]] for x in r])
        split = np.array([x["cases"][-1]["summed_prompt_recoil_with_full_IC"] for x in r])
        defect = float(max(np.abs(C[:,0] - split)))
        ref = max(x["max_time_and_momentum_refinement_difference"] for x in r)
        assert defect < 3e-7 and ref < 2e-6
        return d, C, defect, ref
    def bracket(q, ratio):
        return [[float(q[i]), float(q[i+1])] for i in range(len(q)-1) if (ratio[i]>=.5) and (ratio[i+1]<.5)]
    rows = []
    for f, i in zip(finite_paths, instant_paths, strict=True):
        a, ca, ea, ra = get(f, "ratio5")
        b, cb, eb, rb = get(i, "instant")
        assert a["z"] == b["z"]
        for j, (x,y) in enumerate(zip(a["results"],b["results"],strict=True)):
            assert abs(x["k_h_per_Mpc"] - K[j]) < 1e-12
            for key in ("metric_sha256", "original_native_metric_sha256", "k_h_per_Mpc", "tau_final_Mpc", "initial_k_tau"):
                assert x[key] == y[key], key
        assert len(ca) == len(cb) == len(K)
        assert np.all(ca > 0.) and np.all(cb > 0.) and np.all(ca[:,0] > cb[:,0])
        row = dict(z=a["z"], modes_h_Mpc=K.tolist(), column_meanings=["joined", "prompt", "recoil"],
                   finite=ca.tolist(), instantaneous=cb.tolist(), paired_joined_difference=(ca[:,0]-cb[:,0]).tolist(),
                   finite_downward_half_brackets=bracket(K,ca[:,0]), instant_downward_half_brackets=bracket(K,cb[:,0]),
                   max_joined_minus_split=max(ea,eb), max_refinement=max(ra,rb),
                   shared_metric_sha256=[x["metric_sha256"] for x in a["results"]],
                   input_record_sha256=[sha(f),sha(i)])
        rows.append(row)
    assert sorted(r["z"] for r in rows) == [0.,5.,20.]
    h, C, dh, rh = get(finite_high, "ratio5")
    q = [x["k_h_per_Mpc"] for x in h["results"]]
    assert h["z"] == 20. and np.max(abs(np.array(q)-[26.,30.,34.,40.])) < 1e-12
    return dict(source=source, rows=rows, finite_separate_high_z20=dict(modes_h_Mpc=q, joined=C[:,0].tolist(),
                downward_half_brackets=bracket(q,C[:,0]), max_joined_minus_split=dh, max_refinement=rh,
                input_record_sha256=sha(finite_high)), source_export_sha256=sha(exports),
                recoil_export_sha256=sha(reconstruction),
                scope="Both shapes use identical recorded external gravity and matched redshifts/modes; half intervals are observed downward brackets across sampled modes only. The ratios are scalar density contrasts, not total-matter transfer or a viability test.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    s = p.add_subparsers(dest="mode", required=True)
    r = s.add_parser("recoil")
    for key in ("eos","pairs","directory"):
        r.add_argument(key, type=Path)
    c = s.add_parser("compare")
    c.add_argument("--finite", nargs=3, type=Path, required=True)
    c.add_argument("--instant", nargs=3, type=Path, required=True)
    for key in ("finite_high","reconstruction","exports"):
        c.add_argument(key, type=Path)
    r.add_argument("--sigma", type=float, default=.14)
    a = p.parse_args()
    out = recoil(a.eos,a.pairs,a.directory,a.sigma) if a.mode=="recoil" else compare(a.finite,a.instant,a.finite_high,a.reconstruction,a.exports)
    print(json.dumps(out,indent=2,allow_nan=False))
