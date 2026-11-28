"""Condition the warm pion comparison on Bose and explicit form-factor caps."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import gammaincc,kve
from kloe_pion_band_portal_reference import POINT_MUON, table

from timelike_portal_onsager_coverage_reference import factor,mpi
from three_species_decay_onsager_reference import m,M
from warm_extra_injection_headroom_reference import (
    pion_positive_Q, run_case, VECTOR_MASS, VECTOR_GAMMA,
)


def bose(T):
    return -1./np.expm1(-mpi/T)


def low_strip_per_unit(T,upper_s=.10):
    """Below KLOE, give |Fpi|^2<=1 and the exact p-wave spectral phase."""
    nodes=np.linspace(4*mpi**2,upper_s,185)
    l,h=nodes[:-1],nodes[1:]
    E0,E1=np.sqrt(l),np.sqrt(h)
    r=.25*np.maximum(0.,1-4*mpi**2/h)**1.5
    K0,K1=[kve(j,E0/T)*np.exp(-E0/T) for j in (0,1)]
    W=factor*r*h**1.5*np.diff(nodes)/((h-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2)
    n=q=0.
    for mass,counts in ((m,2),(M,4)):
        z=W*(h+2*mass**2)*np.sqrt(1-4*mass**2/h)
        n+=counts*T*float(z@K1)
        q+=T*float(z@(2*T*K1+E1*K0))
    return np.array([q,n])


def run(eos,measurement,kernel,only=None,return_components=False):
    a=json.loads(Path(kernel).read_bytes())
    K=a["simultaneous_energy_density_coefficient_GeVminus2"]
    TR=.05
    upper_bose=bose(TR)**2
    # Both incoming and final pion may be enhanced in the elastic channel.
    elastic_per_unit_form_factor=upper_bose*gammaincc(5,mpi/TR)/6
    strip_R_max=.25*(1-4*mpi**2/.10)**1.5
    grid,cross,err=table(measurement)
    Rmax=float(max((cross+err)*(grid+.005)/POINT_MUON))
    assert Rmax>0.
    cases=[]
    for label, gap_cap, elastic_cap in (
        ("vacuum_MB_model",0.,0.),
        ("same_band_with_massive_initial_Bose_upper",0.,0.),
        ("Bose_and_threshold_plus_conditional_elastic",1.,1.),
        ("Bose_and_widened_threshold_and_elastic_caps",4.,4.),
    ):
        if only is not None and label != only:
            continue
        is_mb=label.startswith("vacuum")
        def rates(T,data,splits):
            q = pion_positive_Q(T,data,splits)
            if not is_mb:
                # At each fixed T, product of two factors bounded by
                # replacing both incoming energies by the pion rest mass.
                q=(q+gap_cap*low_strip_per_unit(T))*bose(T)**2
            return q
        # The previous coefficient carries Rmax from the measured
        # window. Increase it by the separate bounded-strip maximum.
        tail=1. if is_mb else upper_bose*(1.+gap_cap*strip_R_max/Rmax)
        assert strip_R_max<1.
        z=run_case(eos,K,measurement,980,.008,10,TR=TR,pion_rate_override=rates,
                   pion_tail_multiplier=tail,heat_multiplier=1.+elastic_cap*elastic_per_unit_form_factor,
                   return_components=return_components)
        cases.append(dict(policy=label,threshold_form_factor_squared_cap=gap_cap,
                          spacelike_form_factor_squared_cap=elastic_cap,comparison=z))
    return dict(benchmark_max_initial_two_pion_BE_over_MB=float(upper_bose),
                pion_heat_exponent_over_lepton_per_unit_cap=float(elastic_per_unit_form_factor),
                datasets_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,measurement,kernel)],
                scope="Uses the measured band only as an effective two-massive-pion "
                      "vacuum kernel. Does not make the conditional form-factor "
                      "caps into data or constrain other hadronic channels.",
                cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","measurement","kernel"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.measurement,a.kernel),indent=2,allow_nan=False))
