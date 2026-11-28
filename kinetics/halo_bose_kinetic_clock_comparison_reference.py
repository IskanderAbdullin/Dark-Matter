"""Compare local gross, shear, shape and dynamical rates without equating them."""
import argparse
import json
from pathlib import Path
import numpy as np
from halo_bose_orbital_elastic_clock_reference import R0, V0, GYR, background

def run(gross,shear,shape):
    g=json.loads(gross.read_text())['results']
    s=json.loads(shear.read_text())['examples']
    k=json.loads(shape.read_text())['rows']
    out=[]
    for F in (100.,7000.):
        pot,_=background(np.log1p(1/F))
        for b,c in zip(s,k,strict=True):
            assert b['radius_R']==c['radius_R'] and b['central_occupation']==c['central_occupation']
            if F!=b['central_occupation']:
                continue
            r=b['radius_R']
            a=next(q for q in g if q['peak_occupation']==F and q['pericenter_R']==q['apocenter_R']==r)
            om=np.sqrt(pot(r)[1])/r*V0/R0*GYR
            J=a['stimulated_orbit_encounters_ten_Gyr']/10.
            out.append(dict(central_cap=F,radius_R=r,circular_frequency_per_Gyr=om,
                            true_gross_at_circular_velocity_over_circular_frequency=J/om,
                            local_shear_rayleigh_over_circular_frequency=b['shear_Rayleigh_per_Gyr']/om,
                            local_shape_rayleigh_over_circular_frequency=c['shape_Rayleigh_per_Gyr']/om))
    return dict(scope="Three inequivalent local diagnostics at the same fixed background and cross section. "
                      "Circular frequency is Vc/r in radians per unit time, not the inverse orbital period. "
                      "These are not global modes or a hydrodynamic closure test by themselves.",rows=out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('gross','shear','shape'):
        p.add_argument(a,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.gross,a.shear,a.shape),indent=2,allow_nan=False))
