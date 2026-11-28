"""Refined warm fermion lower bounds with an auditable charged-lepton loss allowance."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from warm_late_visible_recycling_envelope_reference import COUP2,VECTOR_MASS


def read(p):
    return json.loads(Path(p).read_bytes())


def subpole(energy):
    # Sum three final lepton species, with beta_l(1+2m_l^2/s)<=1
    # and (1+2m_dark^2/s)<=1.5. Omit final Pauli and propagator width.
    s=energy**2
    return 3*COUP2*s/(4*np.pi*(VECTOR_MASS**2-s)**2)


def run(kernel,prior,visible):
    K=read(kernel)["max_pair_K_GeVminus2"]
    P,V=read(prior),read(visible)
    goal=P["target"]
    rows=[]
    for p,v in zip(P["cases"],V["cases"],strict=True):
        assert p["q"]==v["q"] and abs(p["anchor_T_GeV"]/v["anchor_T_GeV"]-1.)<1e-12
        n,omega=p["core_start_fermion_number_GeV3"],p["conservative_initial_fermion_Omega"]
        clock=p["integration_time_energy_assisted"]
        E=p["known_energy_upper_GeV4"]+v["additional_early_rho_budget_GeV4"]
        cuts=np.linspace(2,9.9,80)
        cutK=subpole(cuts)
        high_credits=2*E/(n*cuts)
        data=[]
        def shadow(F,k_vis=0.):
            return 1/(1+.5*n*clock*(K*(1+F)**2+k_vis))
        def mixed(F):
            options=shadow(F,cutK)-high_credits
            i=np.argmax(options)
            return max(0.,omega*options[i]),i
        for F in (0.,.1,.5,1.,2.,3.,4.,5.,6.,10.,18.):
            lo,i=mixed(F)
            isolate=float(omega*shadow(F))
            assert 0<=lo<=isolate
            data.append(dict(F_upper=F,isolated_refined_Omega_lower=isolate,
                             with_bidirectional_vacuum_leptons_Omega_lower=float(lo),
                             chosen_energy_cut_GeV=float(cuts[i]),
                             high_energy_budget_as_number_fraction=float(high_credits[i]),
                             low_energy_vector_sigma_v_upper=float(cutK[i])))
        threshold=float(np.sqrt((omega/goal-1)/(.5*n*clock*K))-1)
        mix=brentq(lambda F:mixed(F)[0]-goal,0.,60.,xtol=1e-12)
        assert 1.<mix<threshold
        rows.append(dict(q=p["q"],available_energy_upper_GeV4=float(E),
                         uniform_F_cap_min_to_permit_012_isolated=threshold,
                         uniform_F_cap_min_to_permit_012_vacuum_leptons=mix,
                         data=data))
    return dict(scope="Specified warm network, conditional outgoing scalar occupation cap. "
                      "Allows all known exterior, arbitrary redistribution, "
                      "charged-lepton pair injection, elastic heating and reverse vacuum "
                      "pair loss; no other visible/hadronic loss, new channels or new clocks.",
                target=goal,coefficient=K,cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("kernel","prior","visible"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.kernel,a.prior,a.visible),indent=2,allow_nan=False))
