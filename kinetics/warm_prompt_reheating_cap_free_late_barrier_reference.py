"""Conditional late stable bound without any pointwise scalar occupancy cap."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from three_species_decay_onsager_reference import m, mp, M, y
from warm_prompt_reheating_cold_band_loss_reference import vector_energy_ceiling


def controls():
    rr=[]
    for p,k in ((.82,.16),(.16,.82),(1.21,3.42),(.95,.95)):
        exact=1/max(p,k)
        calculated=quad(lambda x: .5/np.sqrt((p-k)**2+2*p*k*(1-x)),
                        -1.,1.,epsabs=1e-11)[0]
        rr.append(abs(exact/calculated-1))
    ss=np.geomspace(4*m*m,1e6,8000)
    inv=np.sqrt(1-4*m*m/ss)/ss
    peak=1/(6*np.sqrt(3)*m*m)
    assert np.max(inv)<peak
    assert max(rr)<3e-11
    return dict(max_shell_angular_identity_defect=max(rr),
                sampled_energy_weighted_inverse_over_analytic=float(np.max(inv)/peak))


def run(floor,energy,budget):
    f,e,b=[json.loads(p.read_text())['anchors'] for p in (floor,energy,budget)]
    k0=json.loads(floor.read_text())['invariant_vacuum_upper_GeVminus2']
    ka=4*y**4/(M*M)
    kb=ka/(3*np.sqrt(3)*np.pi*m*m)
    cp=(3*np.pi**2)**(2/3)/(2*np.pi**2)
    kv,_=vector_energy_ceiling()
    rows=[]
    for ds,de,db in zip(f,e,b,strict=True):
        T=ds['T_visible_GeV']
        assert T == de['visible_T'] == db['T_visible_GeV']
        H=de['photon_H_floor_at_anchor']
        rho=db['remaining_energy_ceiling']
        A=ds['exposure_coefficient_per_Bose_multiplier']
        n0=2*H*A/k0
        r=rho/(m*n0)
        B=min(2*kv*rho/H,ds['uniform_vector_exposure_coefficient'])
        D=kb*(m*n0)**2/(H*n0)
        for E in (.0003,.001,.002,.004,.0049):
            L=np.sqrt(1-4*E*E/m**2)
            emin=m*m/(4*E)
            r2=1-(m/emin)**2
            dsft=np.pi*ka*(m*n0/mp)/(2*r2*L*emin**3)*n0/(2*H)
            dhrd=np.pi*ka*cp**2/(2*4**(2/3))*(m*n0/E**2)*n0**(1/3)/H
            # Assume X in [cut, 1] before a first downward crossing.
            # Grant at every X the largest scalar energy and light count
            # separately over that interval; this avoids solving a stiff ODE.
            cut=.8
            R=r-cut
            assert R>0
            V=A+B+(dsft+dhrd)*R+D*R*R
            bound=1-V if V<1-cut else None
            rows.append(dict(visible_T_GeV=T,spectral_split_GeV=E,
                             total_energy_to_initial_stable_rest_energy=r,
                             vacuum_quadratic=A,leptonic_loss_upper_on_interval=B,
                             linear_soft_per_residual_energy=dsft,
                             linear_hard_per_residual_energy=dhrd,
                             quadratic_scalar_per_residual_energy_squared=D,
                             comparison_cut=cut,max_loss_on_cut_to_one_interval=V,
                             retained_fraction_lower=bound,
                             omega_h2_lower=bound*ds['initial_Omega_h2_if_counts_freeze']
                             if bound is not None else None))
    assert all(x['retained_fraction_lower'] is None for x in rows[:5])
    assert rows[-1]['omega_h2_lower']>138
    return dict(scope="Anchor-conditional homogeneous factorized isotropic resolved Boltzmann "
                      "collision model for stable pairs into scalars and specified vacuum leptons. "
                      "No scalar occupancy cap; uses Pauli f<=1 for each stable spin. "
                      "No coherent field, novel energy input or additional stable loss channels.",
                invariant_checks=controls(),anchors=rows)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('floor','energy','budget'):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.floor,a.energy,a.budget),indent=2,allow_nan=False))
