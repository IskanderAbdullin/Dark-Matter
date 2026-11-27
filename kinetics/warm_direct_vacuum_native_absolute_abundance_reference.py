"""Use the absolute direct scalar yield instead of diluting the held dark composition."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy import constants
from archived_collisionless_streaming_reference import Background
from conditional_relic_normalization_reference import conversion
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import digest,read
from warm_prompt_reentry_shape_audit_reference import get_mom


def prepare(eos,archive,envelope,higher,out):
    D=read(envelope)
    d=read(archive);S=d["snapshot"]
    assert d["snapshot_sha256"]==digest(S)
    bg=Background(eos)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    T=S["T_visible_GeV"]
    s=(2*np.pi**2/45)*np.exp(bg.gs(np.log(T)))*T**3
    N=get_mom(p,[S["entropy_variable_"+k] for k in ("scalar","fermion")])[0].sum(axis=1)
    Y=N/s
    coef=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    # Independent SI check of the comoving yield normalization.
    s_today_T3=(2*np.pi**2/45)*bg.gs0
    s_today_si=s_today_T3*(bg.T0*1e9*constants.e/(constants.hbar*constants.c))**3
    rho_per_h2=3*(100000/(1e6*constants.parsec))**2*constants.c**2/(8*np.pi*constants.G)
    independent=(1e9*constants.e)*s_today_si/rho_per_h2
    assert abs(independent/coef-1)<1e-4
    held_dims=Y*np.array([mp,m])*coef
    direct=float(D["direct_numerical_Omega_h2"])
    truth=direct/float(direct+held_dims[1])
    dil=direct/(.12*truth)
    assert .001<direct/.12<.002 and 2e-6<truth<2.2e-6
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for tier in ("compact","dense"):
        original=(higher/(tier+"-physical.ini")).read_text()
        for tag,label in ((direct,"absolute"),(D["direct_known_plus_unit_high_Omega_h2_upper"],"upper_normalization")):
            if tier=="dense" and label!="absolute":continue
            q={}
            for l in original.splitlines():
                x,y=l.split("=",1);q[x.strip()]=y.strip()
            q["omega_ncdm"]=f"{tag:.17g}"
            q["omega_cdm"]=f"{(.12-tag):.17g}"
            k=tier+"-"+label
            q["root"]=str(out.resolve()/(k+"_"))
            (out/(k+".ini")).write_text("".join(f"{x} = {y}\n" for x,y in q.items()))
            rows.append(dict(name=k,omega_ncdm=tag,omega_cdm=.12-tag,scalar_fraction_of_0p12=tag/.12))
    return dict(scope="Different normalization scenarios. Direct absolute yield with prescribed stable CDM fill "
                     "keeps scalar abundance derived from only the selected vacuum channels. The old "
                     "held-count ratio instead represents common dilution to imposed dark density. "
                     "Neither gives a fully coupled production model.",
                independent_SI_relative_conversion_difference=float(independent/coef-1),omega_h2_of_held_scalar_and_stable=list(held_dims),held_species_yields=Y.tolist(),
                direct_absolute_omega_h2=direct,direct_plus_held_scalar_fraction=truth,
                dilution_factoring_old_fraction_to_reference_total=float(dil),native_cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","archive","envelope","higher","out"):
        p.add_argument(k,type=Path)
    a=p.parse_args();print(json.dumps(prepare(a.eos,a.archive,a.envelope,a.higher,a.out),indent=2,allow_nan=False))
