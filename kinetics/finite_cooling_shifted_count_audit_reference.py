"""Reevaluate the forward count at the old fixed-loss shifted temperatures."""
import argparse
import json
from pathlib import Path

import numpy as np

from conditional_relic_normalization_reference import conversion
from finite_cooling_frozen_count_retarget_reference import sha
from finite_cooling_reheat_portal_reference import SourceFold, TF
from three_species_decay_onsager_reference import m, mp


def read(path):
    return json.loads(path.read_text())


def prepare(eos,retarget):
    E=sha(eos)
    data=read(retarget)
    assert E in data["input_sha256"].values()
    K=SourceFold(eos,(10,26,20),230,6)
    f=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    rows=[]
    for r in data["controls"][-1]["rows"]:
        primary=r["models"][0]
        assert primary["case"]=="independent_forward"
        t=primary["shifted_root_GeV"]
        ratio=r["T_max_over_T_R"]
        A,B=K.post(t),K.before(t,ratio)
        V=A+B
        gross=f*m*V[1]
        assert abs(gross-.120-primary["stable_omega_depletion"])<2e-11
        rows.append(dict(T_R_GeV=t,Tmax_GeV=t*ratio,scalar_prompt_fixed_count_Omega_h2=f*mp*V[0],
                         expected_gross_stable_Omega_h2=gross,combined_rms_momentum_over_1MeV=(np.sqrt(V[2:]/V[:2])/TF).tolist(),
                         pre_fraction_of_scalar_and_stable_number=(B[:2]/V[:2]).tolist()))
    return dict(controls=[dict(rows=rows)],eos_sha256=E,target_stable_Omega_h2=.120)


def compare(retarget,manifest,old_early,new_early,new_dist,new_late):
    paths=[retarget,manifest,old_early,new_early,new_dist,new_late]
    R,M,E,F,D,L=map(read,paths)
    assert M["eos_sha256"]==E["eos_sha256"]==F["eos_sha256"]==D["eos_sha256"]==L["eos_sha256"]
    assert F["reference_sha256"]==D["reference_sha256"]==sha(manifest) and L["archive_sha256"]==sha(new_dist)
    assert sha(old_early) in R["input_sha256"].values() and D["anchor_temperature_GeV"]==TF
    o=[]
    for r,k,p,q,d,l in zip(R["controls"][-1]["rows"],M["controls"][-1]["rows"],E["episodes"],F["episodes"],
                          D["cases"],L["cases"],strict=True):
        x=k["T_R_GeV"]
        assert x==q["T_R_GeV"]==d["T_R_GeV"]==l["T_R_GeV"]==r["models"][0]["shifted_root_GeV"]
        assert p["T_R_GeV"]==r["archived_unconverted_root_GeV"]
        for v in (q,d,l):
            assert v["Tmax_GeV"]==k["Tmax_GeV"]
        def forward(a):
            return float(np.sum(np.asarray(a["new_scalar_omega"])[:,0]))
        changes=[]
        for a,b in zip(p["modes"],q["modes"],strict=True):
            assert (a["steps"],a["incoming_momentum_bins"])==(b["steps"],b["incoming_momentum_bins"])
            changes.append(forward(b)/forward(a)-1.)
        assert max(changes)-min(changes)<5e-7
        before,after=map(forward,(p["modes"][-1],q["modes"][-1]))
        gross=d["species"][1]["Omega_h2"]
        gross_ref=k["expected_gross_stable_Omega_h2"]
        assert abs(gross/gross_ref-1.)<1e-6
        late=l["controls"][-1]["perturbative_unblocked_scalar_Omega_h2_equivalent"]
        # This normalization control preserves the new spectrum shape and
        # clock. The late forward count scales quadratically in density.
        early_loss=m/mp*after
        corrected_late=late*(1-early_loss/gross)**2
        net=gross_ref-m/mp*(after+corrected_late)
        omega_residual=net-M["target_stable_Omega_h2"]
        dT=-omega_residual/r["stable_count_temperature_slope_per_GeV"]
        assert abs(dT)*1e6<.05
        o.append(dict(ratio=r["T_max_over_T_R"],test_root_GeV=x,prompt_stable_omega_archive_defect=gross-gross_ref,
             prior_gross_early_scalar_omega=before,new_gross_early_scalar_omega=after,
             relative_change_early_forward_by_matched_quadrature=changes,new_gross_late_scalar_omega=late,
             previous_and_new_scalar_rms_ptf=[v["modes"][-1]["first_order_mixture_scalar_rms_p_over_TF"] for v in (p,q)],
             late_after_density_only_control_scalar_omega=corrected_late,
             recovered_stable_omega=net,residual_stable_omega=omega_residual,
             one_step_temperature_update_keV=dT*1e6))
    return dict(cases=o,record_sha256={p.name:sha(p) for p in paths},
                scope="Recomputed first-order forward conversion at the fixed-loss displaced roots. Prompt stable spectra and the radiation clock are regenerated; optional late density-only attenuation does not change the distribution shape. The final one-step temperature is a slope proxy, not a jointly iterated signed kinetic solution.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    s=p.add_subparsers(dest="mode",required=True)
    a=s.add_parser("prepare")
    for k in ("eos","retarget"):
        a.add_argument(k,type=Path)
    b=s.add_parser("compare")
    for k in ("retarget","manifest","old_early","new_early","new_dist","new_late"):
        b.add_argument(k,type=Path)
    a=p.parse_args()
    out=prepare(a.eos,a.retarget) if a.mode=="prepare" else compare(a.retarget,a.manifest,a.old_early,a.new_early,a.new_dist,a.new_late)
    print(json.dumps(out,indent=2,allow_nan=False))
