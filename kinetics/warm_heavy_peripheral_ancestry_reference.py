"""Kinematic support envelope for a finite number of retained light inputs."""

import argparse
import json
import math
from pathlib import Path


from three_species_decay_onsager_reference import M, m, mp


def envelope(P, leaves, b):
    """Maximum squared comoving momentum allowed by the energy budget."""
    return leaves**2*P**2+((leaves*m)**2-M**2)*b**2


def verify_redshift_step(P):
    # If a light/scalar intermediate was made from k core leaves at early b,
    # its momentum stays fixed until later b. Its effective energy allowance
    # redshifts consistently because out_mass <= k*m.
    for k in (1,2,3,4,5):
        for out_mass in (m,mp):
            for b0,b1 in ((1.,1.),(1.,1.02),(1.,100.),(15.,1000.)):
                p2 = k*k*P*P+((k*m)**2-out_mass**2)*b0*b0
                after = p2+(out_mass*b1)**2
                allowance = (k*math.hypot(P,m*b1))**2
                assert after <= allowance*(1+3e-15)


def run(path):
    data=json.loads(path.read_text())
    rows=[]
    for c in data["cases"]:
        P=c["assumed_equal_comoving_daughter_cap_float_GeV"]
        cut=c["high_cut_GeV"]
        verify_redshift_step(P)
        limits=[]
        for leaves in range(1,9):
            C=(leaves*m)**2-M**2
            bound=math.sqrt(max(0.,envelope(P,leaves,1.)))
            # For C<=0 the largest value is the initial one. For C>0 and
            # bound<cut, record the latest b before the sufficient exclusion
            # can fail. Lack of exclusion is not kinematic attainability.
            b_cross=math.sqrt((cut**2-leaves**2*P**2)/C) if C>0 and bound<cut else None
            limits.append(dict(core_leaves=leaves, initial_heavy_ceiling_GeV=bound,
                               exclude_all_b_ge_1=(C<=0 and bound<cut),
                               b_cross_if_initially_excluded=b_cross))
        assert limits[2]["exclude_all_b_ge_1"]
        if c["q"]==16:
            assert limits[4]["exclude_all_b_ge_1"] and not limits[5]["exclude_all_b_ge_1"]
        if c["q"]==32:
            assert limits[3]["initial_heavy_ceiling_GeV"]>cut
        rows.append(dict(q=c["q"],P_GeV=P,high_cut_GeV=cut,controls=limits))
    return dict(scope="All core leaves light or scalar on retained support, all propagating "
                      "intermediates light or scalar with only collisionless redshift between "
                      "events, and no energy from other ancestors. Bounds apply to H creation "
                      "only, not subsequent heavy elastic redistribution.",cases=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("support",type=Path)
    print(json.dumps(run(parser.parse_args().support),indent=2,allow_nan=False))
