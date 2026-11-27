"""Use isotropic incoming directions and the scalar count, not a soft Bose cap."""
import argparse
import json
from pathlib import Path

import numpy as np

from archived_scalar_pair_reaction_reference import amplitude
from three_species_decay_onsager_reference import m, mp, M, y
from warm_prompt_reheating_cold_band_loss_reference import retention, vector_energy_ceiling


def angular_controls():
    # Constant |M|^2/s=1 example: exact fixed-magnitude s interval and
    # analytic integral of s/sqrt(E^2-s), evaluated independently by quadrature.
    from scipy.integrate import quad
    defects, ratios = [], []
    for e in (2*mp, 1e-6, 1e-5, 1e-4, 1e-3):
        k = np.sqrt(e*e-mp*mp)
        min_leg = m*m/(4*e)
        for scale in (2.5, 10., 100.):
            e1 = scale*min_leg
            e2 = 1.2*e1
            p1, p2 = np.sqrt(e1*e1-m*m), np.sqrt(e2*e2-m*m)
            E = e1+e2
            slo = 2*m*m+2*(e1*e2-p1*p2)
            shi = 2*m*m+2*(e1*e2+p1*p2)
            # Squared support condition: (s-2Ee)^2 <= 4k^2(E^2-s).
            centre = 2*E*e-2*k*k
            delta = 2*k*np.sqrt(E*E-2*E*e+k*k)
            lo, hi = max(slo,centre-delta), min(shi,centre+delta)
            if hi <= lo:
                continue
            assert hi-lo <= 4*E*k*(1+1e-12)
            assert hi <= 4*E*e*(1+1e-12)
            true = quad(lambda s:s/np.sqrt(E*E-s),lo,hi,epsabs=1e-30)[0]
            Pmin = np.sqrt(E*E-hi)
            Pmax = np.sqrt(E*E-lo)
            # Stable expression for the difference of primitives.
            mixed = (E*E*(lo+hi)-lo*hi)/(E*E+Pmax*Pmin)
            analytic = 2*(hi-lo)/(Pmax+Pmin)*(lo+hi+mixed)/3
            defects.append(abs(analytic/true-1))
            L = np.sqrt(1-4*e*e/(m*m))
            # The stronger previous bound E >= m^2/e holds on contributing intervals.
            assert E >= m*m/e*(1-1e-12)
            ratio = true/(4*e/L*4*E*k)
            assert ratio <= 1+1e-12
            ratios.append(ratio)
    assert len(ratios) >= 8 and max(defects)<1e-11
    s = np.geomspace(4*m*m*(1+1e-8),1e7,200)
    a = amplitude(s[:,None],np.linspace(-1,1,257))
    maxamp = np.max(a/s[:,None])/(4*y**4/(M*M))
    assert maxamp < 1.
    return dict(support_cases=len(ratios), max_quadrature_identity_defect=max(defects),
                max_fraction_of_angular_upper=max(ratios),
                max_sampled_invariant_amplitude_fraction=float(maxamp))


def run(floor, energy, bounds):
    fl = json.loads(floor.read_text())
    en = json.loads(energy.read_text())
    bd = json.loads(bounds.read_text())
    ka = 4*y**4/(M*M)
    k0 = fl['invariant_vacuum_upper_GeVminus2']
    kv, _ = vector_energy_ceiling()
    rows = []
    for d, e, b in zip(fl['anchors'], en['anchors'], bd['anchors'],strict=True):
        T = d['T_visible_GeV']
        assert T == e['visible_T'] == b['T_visible_GeV']
        rho = b['remaining_energy_ceiling']
        H = e['photon_H_floor_at_anchor']
        n = 2*H*d['exposure_coefficient_per_Bose_multiplier']/k0
        n_soft = rho/mp
        bv = 2*kv*rho/H
        v0 = d['uniform_vector_exposure_coefficient']
        raw = d['initial_Omega_h2_if_counts_freeze']
        for fh in (1.,10.):
            for lim in (1e-6,1e-5,1e-4,1e-3):
                emin = m*m/(4*lim)
                assert emin > m
                r2 = 1-(m/emin)**2
                L = np.sqrt(1-4*lim*lim/(m*m))
                extra = (1+fh)*ka*np.pi*n_soft/(2*r2*L*emin**3)
                A = k0*(1+fh)**2*n/(2*H)
                delta = extra*n/(2*H)
                # Grant the largest allowed physical soft density forever.
                # This is weaker than integrating its b^-3 decrease.
                low = raw*max(retention(A+delta,bv),retention(A+delta+v0,0.))
                rows.append(dict(T_visible_GeV=T,assumed_hard_cap=fh,
                                 split_GeV=lim,minimum_energy_per_light_leg_GeV=emin,
                                 conservative_soft_count_ceiling_GeV3=n_soft,
                                 excess_stimulated_sigma_v_upper_GeVminus2=extra,
                                 extra_quadratic_exposure=delta,conditional_light_omega_floor=low))
    return dict(scope="Homogeneous isotropic factorized incoming fermion distributions, "
                      "standard vacuum 2-to-2 Boltzmann collision operator, continuum scalar occupation, "
                      "uniform hard-band cap only, standard photon/entropy comparison and prior "
                      "charged-lepton energy allowance. No constraint on soft occupancy supremum, "
                      "but not a model of a coherent condensate, anisotropic streams, or extra channels.",
                invariant_controls=angular_controls(),rows=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('floor','energy','bounds'):
        p.add_argument(name,type=Path)
    a = p.parse_args()
    print(json.dumps(run(a.floor,a.energy,a.bounds),indent=2,allow_nan=False))
