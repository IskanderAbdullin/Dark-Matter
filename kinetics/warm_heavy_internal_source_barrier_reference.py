"""Use energy moments to bound production at every heavy outgoing momentum."""
import argparse
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import M,m,mp,y,width


def read(p):
    return json.loads(Path(p).read_bytes())


def contributions(E,cap):
    # Uniform overcount of spin, labels and both light initial charge channels.
    label=32.
    phase=label/(16*np.pi**3)  # after p2 delta, p1 angular avg, and tag 1/(2 E_H)
    # Incoming scalars: E1E2>=M^2. If one scalar is soft (E<M/2),
    # its partner is hard (E>2M). The angle integral on the soft side
    # uses 1/max(k,p_H), not 1/k alone.
    beta_h=np.sqrt(1.-(2*mp/M)**2)
    beta_tail=np.sqrt(1.-(mp/(2*M))**2)
    assert (1-(mp/(3*M/8))**2)**(-.5)<2.
    integral=(2*np.pi**2*E)**2/M**4*(16/beta_h**2+4/beta_tail)
    # Amplitude divided by incoming energy denominators <= 2y^4/m^2.
    scalar_job=phase*2*y**4/m**2*integral

    n=E/m
    I=.5*(3*np.pi**2*n)**(2/3)*cap**(1/3)
    J1=I/m
    J2=m*I+np.pi**2*n
    cc=4*((M+m)**2-mp*mp)/(M*M-m*m)**2
    # t,u <= -(M^2-m^2)^2/s. After splitting (1+u)^2,
    # the initial energy measures leave E1E2 and (E1E2)^-1.
    light_job=phase*y**4/2*(J1*J1+(cc*J2)**2)
    return dict(external_E_budget_GeV4=E,
                scalar_pair_upper_in_units_Gamma_M_over_E_parent=float(scalar_job/(width*M)),
                light_pair_upper_in_units_Gamma_M_over_E_parent=float(light_job/(width*M)),
                conservative_spin_and_channel_overcount=label)


def run(pion,visible):
    P,V=read(pion),read(visible)
    rows=[]
    for z in P["cases"]:
        E=z["controls"][-1]["new_available_energy_upper_GeV4"]
        trials=[]
        for c in V["visible_source_barrier_checks"]:
            delta=c["cap"]
            a=contributions(E,delta)
            r=delta*(1-delta)
            q=(a["scalar_pair_upper_in_units_Gamma_M_over_E_parent"]+
               a["light_pair_upper_in_units_Gamma_M_over_E_parent"]+
               r*c["max_visible_source_over_decay_margin"])
            assert 0<q<r
            trials.append(dict(cap=delta,total_bounded_input_over_outgoing_decay_margin=q/r,
                               remaining_source_budget_in_units_Gamma_M_over_E_parent=r-q,
                               **a))
        rows.append(dict(q=z["q"],cases=trials))
    return dict(scope="Modewise heavy-fermion barrier under a cap on the light occupations, "
                      "finite energy and spatial homogeneity/isotropy, with the specified "
                      "number-changing channels and positive budgets. Does not bound gains "
                      "from elastic redistribution of the heavy fermion, or the light cap itself.",
                results=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("pions")
    p.add_argument("visible_barrier")
    a=p.parse_args()
    print(json.dumps(run(a.pions,a.visible_barrier),indent=2,allow_nan=False))
