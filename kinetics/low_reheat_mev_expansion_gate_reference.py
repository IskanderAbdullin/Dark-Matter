"""Extra homogeneous expansion from the represented 1 MeV benchmark."""
import argparse
import json
from pathlib import Path

import numpy as np

from scalar_pair_expanding_continuation_reference import MASS,moments
from late_radiation_background_reference import TabulatedRadiation


def gate(table,early,late,shape):
    x=json.loads(Path(early).read_text())["rows"][-1]
    k=json.loads(Path(late).read_text())["states"][-1]
    s=json.loads(Path(shape).read_text())["snapshots"][-1]
    b=k["scale_factor"]
    ai=x["final_a"]
    T=s["T_visible_GeV"]
    assert abs(T/.001-1.)<2e-6 and abs(s["scale_factor"]/(ai*b)-1.)<2e-12
    p=np.array(s["lab_momentum_edges_GeV"])
    U=moments(p,[s["entropy_variable_scalar"],s["entropy_variable_fermion"]])
    rho=(U[:,1]+MASS[:,None]*U[:,0]).sum(axis=1)
    ref=TabulatedRadiation(table,T)
    r=ref.state(1.)[1]
    assert max(abs(U[:,0].sum(axis=1)*s["scale_factor"]**3/np.array(s["represented_species_count_a3"])-1))<2e-12
    ext_early_at_anchor=x["visible_bath_feedback"]["exterior_a4_energy_distribution_independent_upper_bound"]/ai**4
    ext_mid=k["exterior"]["exterior_a4_energy_distribution_independent_upper_bound"]/b**4
    # Redshift drops each existing particle energy; keeping the energy of
    # every early cohort constant gives a deliberately generous upper.
    re=rho.sum(); up=re+ext_early_at_anchor/b**3+ext_mid
    dh=lambda f:f/(1.+np.sqrt(1.+f))
    nu=np.pi**2*(7./120)*T**4
    return dict(T_photon_GeV=T,late_pair_only_expansion_ratio=b,
        retained_species_energy_relative_to_radiation=(rho/r).tolist(),
        early_exterior_energy_over_radiation_upper=ext_early_at_anchor/b**3/r,
        newly_created_exterior_energy_over_radiation_upper=ext_mid/r,
        fractional_increase_H_including_exteriors_bracket=[dh(re/r),dh(up/r)],
        total_dark_energy_in_units_of_ideal_Tnu_eq_Tgamma_single_flavor=[re/nu,up/nu],
        scope="Only the extra instantaneous homogeneous density from particles in the tested low-reheat benchmark. Prescribed visible bath, no neutrino transport, nuclear network, modified weak rates, or BBN likelihood.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("table","early","late","shape"):p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(gate(a.table,a.early,a.late,a.shape),indent=2,allow_nan=False))
