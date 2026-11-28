"""Compare hypothetical local fast-relaxation requirement with fixed Born quartics."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import constants

CM2_PER_G_PER_GEV_MINUS3=(constants.hbar*constants.c/(1e9*constants.e)*100)**2/(1e9*constants.e/constants.c**2*1000)
M=.000000086
SCALE=1.25


def sigma_per_mass(lam):
    return float(CM2_PER_G_PER_GEV_MINUS3*lam*lam/(128*np.pi*M**3))


def run(kin):
    cases={}
    for lam in (7.65175e-8,1.08515714e-6):
        cases[lam]=sigma_per_mass(lam)
    assert abs(cases[7.65175e-8]-5.)<1e-4
    assert abs(cases[1.08515714e-6]-1005.62)<.03
    rows=[]
    vals=json.loads(kin.read_text())['rows']
    for st in vals:
        the=[]
        for field in ('local_shear_rayleigh_over_circular_frequency',
                      'local_shape_rayleigh_over_circular_frequency'):
            f=st[field]
            need=SCALE/f
            lam=float(np.sqrt(need*128*np.pi*M**3/CM2_PER_G_PER_GEV_MINUS3))
            assert abs(sigma_per_mass(lam)/need-1)<4e-15
            the.append(dict(trial_function=field,
                            necessary_sigma_over_mass_cm2_per_g_for_rate_to_beat_circular_clock=need,
                            necessary_quartic_lambda_for_that_contact_kernel=lam,
                            quotient_over_circular_frequency_at_5_cm2_per_g=f*5/SCALE,
                            quotient_over_circular_frequency_at_static_curvature_lambda=f*cases[1.08515714e-6]/SCALE))
        rows.append(dict(central_cap=st['central_cap'],radius_R=st['radius_R'],data=the))
    return dict(scope="Fixed halo, standard Born real scalar quartic, linear local kinetic rates only. "
                      "Requirements concern rates larger than the specified circular frequency, "
                      "not a galaxy-age equilibrium criterion or a universal observational limit. "
                      "Static plasma curvature benchmark is separate and conditional.",
                Born_cross_section_over_mass={str(k):v for k,v in cases.items()},
                rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('kinetics_clock',type=Path)
    print(json.dumps(run(p.parse_args().kinetics_clock),indent=2,allow_nan=False))
