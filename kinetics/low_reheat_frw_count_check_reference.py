"""Independent full support and underflow checks at the low reheat benchmark."""
import argparse
import json
from pathlib import Path

import numpy as np

from conditional_relic_normalization_reference import conversion, proxy
from empty_lepton_pair_spectrum_reference import Source
from joint_visible_fd_source_reference import Visible
from joint_visible_heavy_adiabatic_reference import shell_projection
from joint_visible_quantum_source_reference import signed_shells
from late_radiation_background_reference import TabulatedRadiation
from residual_visible_pair_production_reference import interval
from sm_dark_entropy_exchange_reference import StandardModel
from three_species_decay_onsager_reference import M,m
from visible_frw_energy_ledger_reference import ConservativeEOS


def cold_reaction_regression():
    T=.001127
    edges=np.r_[0.,np.geomspace(.002,5.,20)]
    order=(6,12,14)
    z=signed_shells(T,lambda e:np.full_like(e,np.inf),edges)
    d=signed_shells(T,lambda e:e/T,edges)
    raw=Visible(order,polar=16).rates(T,edges,1.)
    direct=np.sum(abs(z[0]-raw[0]),axis=1)/raw[0].sum(axis=1)
    balance=np.sum(abs(d[0]-d[1]),axis=1)/raw[0].sum(axis=1)
    assert max(direct)<3e-12 and max(balance)<3e-13
    assert z[1].sum()==0.
    heavy=shell_projection(T,np.zeros_like,np.zeros_like,edges,order=order)
    dc=np.sum(abs(heavy[:,0]-raw[[2,1]]),axis=2)/raw[[2,1]].sum(axis=2)
    assert dc.max()<3e-12 and heavy[:,1].sum()==0.
    # Only q>2 GeV carries this inverse current. There the visible
    # exp(-q/T) is exactly zero in float64; the inverse must survive.
    hot=signed_shells(T,lambda e:np.where(e>1.,np.log(3.),np.inf),edges)
    assert hot[1,0].sum()>1e-16
    assert np.max(abs(hot[0]-z[0]))==0.
    Qmax=Source(order,dark_mass=m).geometry(T)[1].max()/T
    assert Qmax>750
    return dict(T_GeV=T,maximum_node_energy_over_temperature=Qmax,
                relative_cold_vacuum_direct_L1=direct.tolist(),
                relative_cold_vacuum_heavy_L1=dc.tolist(),
                relative_cold_KMS_direct_L1=balance.tolist(),
                manufactured_direct_inverse_count_both_incoming_energies_above_1_GeV=hot[1,0].sum())


def compare(table,files):
    k=conversion(table)["Omega_h2_per_GeV_times_yield"]
    eos=ConservativeEOS(StandardModel(table,lower_knot_limit=.001))
    out=[]
    for name in files:
        for r in json.loads(Path(name).read_text())["rows"]:
            tr=r["initial_visible_T_GeV"]; tf=r["final_visible_temperature_GeV"]
            clock=TabulatedRadiation(table,tr)
            sf=eos.state(tr)[2]+r["visible_bath_feedback"]["visible_comoving_entropy_change"]
            kept=np.array(r["final_counts_div_initial_visible_entropy"])*eos.state(tr)[2]/sf
            outside=(np.array(r["comoving_source_channel_count_and_total_energy_tails"])
                     +r["first_pulse_count_and_total_energy_tails"])[:,0].sum(axis=1)/sf
            allp=kept+np.array([outside[2],outside[:2].sum()])
            direct,heavy=interval(clock,tf,tr,(11,28,18),6)
            dleft,hleft=interval(clock,.001,tf,(11,28,18),6)
            prescribed=np.array([heavy,direct+heavy])
            endpoints=np.geomspace(tf,tr,101)
            relrho=max(abs(eos.state(t)[0]/(np.pi**2/30*t**4*np.exp(clock.rg(np.log(t))))-1.)
                       for t in endpoints)
            assert relrho<2e-3
            defect=allp/prescribed-1.
            assert abs(defect[1])<.04
            out.append(dict(T_R_GeV=tr,T_final_GeV=tf,steps=r["nsteps"],
                state_dependent_source=r["state_dependent_quantum_visible_source"],
                retained_fixed_count_Omega_h2=proxy(kept,k),
                retained_plus_source_tails_fixed_count_Omega_h2=proxy(allp,k),
                independent_prescribed_full_vacuum_Omega_h2=proxy(prescribed,k),
                prescribed_remaining_vacuum_Omega_h2_to_1_MeV=proxy([hleft,dleft+hleft],k),
                relative_complete_minus_prescribed_vacuum=defect.tolist(),
                relative_outside_yield_to_all=(1-kept/allp).tolist(),
                maximum_relative_conservative_vs_tabulated_radiation_density=relrho,
                relative_visible_change_at_same_a=r["visible_bath_feedback"]["relative_visible_temperature_shift"]))
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("table"); p.add_argument("frw",nargs="+")
    a=p.parse_args()
    print(json.dumps(dict(cold_regression=cold_reaction_regression(),controls=compare(a.table,a.frw)),
                     indent=2,allow_nan=False))
