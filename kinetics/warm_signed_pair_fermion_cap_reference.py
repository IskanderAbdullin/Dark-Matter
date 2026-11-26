"""Conditional fermion survival without a scalar occupation cap."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from three_species_decay_onsager_reference import m,M,y
from warm_fermion_survival_with_leptons_reference import subpole


def read(p):
    return json.loads(Path(p).read_bytes())


def run(fermion,leptons,isolated):
    F,L,I=map(read,(fermion,leptons,isolated))
    K=L["coefficient"]; goal=F["target"]
    m_i=np.array([m,M])
    m_x=np.array([M,m])
    etas=np.geomspace(1.005,6.,260)
    S=np.max(2*y**4/(np.pi*m_x**2*m_i**4)) / np.sqrt(1-etas**-2)
    H=np.max(2*y**4/(np.pi**3*m_x**2*m_i**2))*etas**2
    assert min(S)>0 and min(H)>0
    cuts=np.linspace(2.,9.9,80)
    U=subpole(cuts)
    all=[]
    for a,b,c in zip(F["cases"],L["cases"],I["cases"],strict=True):
        q=a["q"]
        assert b["q"]==c["q"]==q
        N=a["core_start_fermion_number_GeV3"]; O=a["conservative_initial_fermion_Omega"]
        E=b["available_energy_upper_GeV4"]
        I3=a["integration_time_energy_assisted"]
        clock=c["bounds"][1]["clock"]
        stop=clock["b_stop"]; h=clock["floor_min"]; Hp=clock["photon_anchor_GeV"]
        def moment(p):
            return ((1-stop**(2-p))/h + stop**(2-p))/(Hp*(p-2))
        I4,I5=moment(4),moment(5)
        # Nonstimulable reference quadratic; both quadratic rates use
        # n_plus*n_minus<=n_species^2/4. Above-cut lepton losses are reserved.
        survivors=1/(1+.5*N*I3*(K+U))-2*E/(N*cuts)
        j=np.argmax(survivors)
        def conditional(f):
            # Fermi-ball rearrangement per charge: I_rad=∫p f dp.
            Jmax=.5*(3*np.pi**2*N)**(2/3)*f**(1/3)
            Jhalf=Jmax*2**(-2/3)
            # At the final descent below the retained initial inventory,
            # F_fermion <= N. The omitted stimulation is a positive debit.
            dsoft=2*S*Jmax*E*E/4*I5 / N
            dhard=2*H*Jhalf**2*E*I4 / N
            k=np.argmin(dsoft+dhard)
            z=float(O*max(0,survivors[j]-dsoft[k]-dhard[k]))
            return dict(fermion_per_spin_cap=f,scalar_per_mode_cap=None,
                        lower_fermion_Omega_h2=z,eta=float(etas[k]),
                        additive_soft_and_hard_loss_fraction=[float(dsoft[k]),float(dhard[k])],
                        lepton_split_GeV=float(cuts[j]))
        rows=[conditional(f) for f in (0.,.0001,.001,.003,.01,.02,.05,.1,.5)]
        assert abs(rows[0]["lower_fermion_Omega_h2"]/b["data"][0]["with_bidirectional_vacuum_leptons_Omega_lower"]-1.)<2e-13
        threshold=brentq(lambda f:conditional(f)["lower_fermion_Omega_h2"]-goal,0.,.5,xtol=2e-15)
        all.append(dict(q=q,additional_moment_I4_I5=[I4,I5],finite_fermion_cap_needed_to_evade=threshold,
                        controls=rows))
    return dict(scope="Spatially homogeneous isotropic vacuum signed dark pair network, "
                      "charged-lepton sources, losses and heat as earlier. Uses per-spin cap "
                      "on light and heavy fermion occupations, including any reentered tails. "
                      "No bound on scalar occupancy, no claim of a dynamical fermion cap.",
                target=goal,cases=all)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for t in ("fermion","leptons","isolated"):
        p.add_argument(t)
    a=p.parse_args()
    print(json.dumps(run(a.fermion,a.leptons,a.isolated),indent=2,allow_nan=False))
