"""Distinguish true Yukawa on-shell support from sparse Gauss event misses.

Uses exactly the quadratic scalar energy-cell convention and 12x12
Gauss parent-angle events of phi_bose_cell_entropy.fixture.
"""
import math
from numpy.polynomial.legendre import leggauss
from kinetics.phi_bose_ir_accessibility import (
    M1, M2, MPHI, rest_kinematics, min_over_parent_cutoff,
)


def sampled_minimum(pmax, nr=12, nangle=12):
    x,_=leggauss(nr)
    c,_=leggauss(nangle)
    estar,qstar=rest_kinematics()
    values=[]
    for xx in x:
        p=pmax*(xx+1)/2
        e1=math.hypot(p,M1)
        for cc in c:
            values.append((e1*estar+p*qstar*cc)/M1)
    return min(values)


def first_cell_edge(pmax,n):
    if pmax<=0 or n<2: raise ValueError("Positive cutoff and >=2 cells")
    return MPHI+(math.hypot(pmax,M1)+.2-MPHI)/n**2


def audit():
    out=[]
    for pmax in (2.,8.,16.):
        true_min=min_over_parent_cutoff(pmax)
        sample_min=sampled_minimum(pmax)
        for n in (8,16,32,64):
            edge=first_cell_edge(pmax,n)
            out.append({"pmax_GeV":pmax,"n":n,"physical_min_GeV":true_min,
                        "sampled_min_GeV":sample_min,"first_cell_top_GeV":edge,
                        "physical_support_in_first_cell":bool(edge>true_min),
                        "sampled_event_in_first_cell":bool(edge>sample_min),
                        "quadrature_false_zero":bool(true_min<edge<=sample_min)})
    return out


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2))
