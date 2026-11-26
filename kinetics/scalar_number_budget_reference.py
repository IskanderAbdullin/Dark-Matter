"""Scalar rest-mass upper bound conditional on late total-count conservation."""
import argparse
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m, mp


def budget(name,omega,reference):
    # Both opposite fermion charges are included in the stable Omega input.
    phi,chi=map(float,omega)
    r=mp/m
    all_scalar=phi+r*chi
    assert phi<=all_scalar<=phi+chi and all_scalar>0.
    x=(phi+chi-reference)/(1-r)
    final=None
    if 0.<=x<=chi:
        final=[phi+r*x,chi-x]
        assert abs(sum(final)/reference-1.)<2e-13
        assert abs((final[0]+r*final[1])/all_scalar-1.)<2e-15
    return dict(case=name,initial_fixed_count_Omega_h2=[phi,chi],
        scalar_upper_Omega_h2_if_all_particles_become_scalars=all_scalar,
        scalar_ceiling_fraction_of_reference=all_scalar/reference,
        scalar_number_multiplier_needed_for_half_reference=.5*reference/all_scalar,
        scalar_number_multiplier_needed_for_99_percent_reference=.99*reference/all_scalar,
        scalar_number_multiplier_needed_for_all_reference=reference/all_scalar,
        conditional_species_Omega_h2_if_conversion_only_matches_total_reference=final)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("normalization")
    a=p.parse_args()
    d=json.loads(Path(a.normalization).read_text())
    ref=d["Planck_cold_dark_matter_Omega_h2_reference"]
    rows=[budget(r["archive"],r["conditional_Omega_h2"],ref) for r in d["at_fixed_archive_counts"]]
    sweep=d["prescribed_vacuum_sweep"]
    r=next(r for r in sweep["rows"] if r["T_R_GeV"]==.05)
    rows.append(budget("50 MeV full prescribed vacuum source to 1 MeV",r["conditional_Omega_h2"],ref))
    rows.append(budget("9.019 MeV prescribed vacuum root",sweep["fine_proxy_at_nominal_T"],ref))
    assert np.max([r["scalar_ceiling_fraction_of_reference"] for r in rows])<.005
    print(json.dumps(dict(m_scalar_over_m_fermion=mp/m,observational_reference_Omega_h2=ref,
         scope="Fixed initial total dark number per visible comoving entropy; subsequent conversion among these retained species can preserve or reduce the total, not increase it. No additional sources, number multiplication, cooling and clustering assumptions, or applicability to other model parameters.",
         rows=rows),indent=2,allow_nan=False))
