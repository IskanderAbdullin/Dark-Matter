"""Number-conserving relic proxies, not a prediction with depletion."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy import constants
from scipy.optimize import brentq

from late_radiation_background_reference import TabulatedRadiation
from residual_visible_pair_production_reference import interval
from sm_dark_entropy_exchange_reference import StandardModel
from three_species_decay_onsager_reference import m,mp


def conversion(table,T0=2.72548):
    g=float(np.loadtxt(table)[0,1])
    entropy=2*np.pi**2/45*g*(constants.k*T0/(constants.hbar*constants.c))**3
    H100=100000./(1e6*constants.parsec)
    critical=3*H100**2/(8*np.pi*constants.G)
    unit=1e9*constants.e/constants.c**2
    return dict(CMB_temperature_K=T0,late_standard_model_entropy_degrees=g,
                entropy_cm3=entropy/1e6,critical_kg_per_m3_over_h2=critical,
                Omega_h2_per_GeV_times_yield=entropy*unit/critical)


def proxy(y,factor):
    return (factor*np.array([mp,m])*np.asarray(y)).tolist()


def archives(eos,rows,factor,reference=.120):
    out=[]
    for label,r in rows:
        y=np.array(r["final_counts_div_initial_visible_entropy"])
        s0=eos.state(r["initial_visible_T_GeV"])[2]
        visible=r["visible_bath_feedback"]
        entropy_factor=1+visible["visible_comoving_entropy_change"]/s0
        y=y/entropy_factor
        val=proxy(y,factor)
        out.append(dict(archive=label,final_visible_temperature_GeV=r["final_visible_temperature_GeV"],
                        yield_per_retained_species=y.tolist(),visible_entropy_ratio=entropy_factor,
                        conditional_Omega_h2=val,stable_proxy_over_Planck=val[1]/reference,
                        stable_depletion_fraction_if_no_extra_dilution=max(0.,1.-reference/val[1])))
    return out


def sweep(table,factor):
    bg=TabulatedRadiation(table,.020)
    def sample(T,order=(8,18,14),time=5):
        d,h=interval(bg,.001,T,order,time)
        return proxy((h,d+h),factor)
    samples=[]
    for T in (.003,.005,.007,.009,.012,.015,.018,.020,.024,.03,.05):
        o=sample(T)
        samples.append(dict(T_R_GeV=T,conditional_Omega_h2=o))
    assert np.all(np.diff([q["conditional_Omega_h2"][1] for q in samples])>0)
    a,b=next((a,b) for a,b in zip(samples,samples[1:])
             if a["conditional_Omega_h2"][1]<.120<b["conditional_Omega_h2"][1])
    T=brentq(lambda u:sample(u)[1]-.120,a["T_R_GeV"],b["T_R_GeV"],xtol=3e-9)
    fine=sample(T,(11,28,18),8)
    assert abs(fine[1]/.120-1.)<3e-5
    return dict(rows=samples,vacuum_only_nominal_T_reheat_for_Planck_stable_count_GeV=T,
                fine_proxy_at_nominal_T=fine,
                scope="prescribed radiation, vacuum direct+prompt pair production only; no inverse/depletion or visible feedback")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos");p.add_argument("radiation");p.add_argument("archives",nargs="+")
    p.add_argument("--no-sweep",action="store_true")
    a=p.parse_args()
    c=conversion(a.radiation)
    f=c["Omega_h2_per_GeV_times_yield"]
    files=[(Path(x).name,json.loads(Path(x).read_text())["rows"][-1]) for x in a.archives]
    result=dict(conversion=c,Planck_cold_dark_matter_Omega_h2_reference=.120,
                at_fixed_archive_counts=archives(StandardModel(a.eos),files,f))
    if not a.no_sweep:
        result["prescribed_vacuum_sweep"]=sweep(a.radiation,f)
    print(json.dumps(result,indent=2,allow_nan=False))
