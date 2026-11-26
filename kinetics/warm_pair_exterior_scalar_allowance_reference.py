"""Sufficient bound on additional exterior Bose occupancy for the frozen pair sign."""

import argparse
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import M,m,mp,y
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_frozen_continuum_cap_bound_reference import bound


def allowance(path):
    info=bound(path)
    a=read(path)
    a=a[-1] if isinstance(a,list) else a
    S=a["snapshot"]
    assert digest(S)==info["snapshot_sha256"]
    d=np.asarray(S["entropy_variable_scalar"])
    zmin=float(min(np.min(d[0]),np.min(d[0]+d[1])))
    assert zmin>0.
    bmax=float(np.exp(-zmin)/(-np.expm1(-zmin)))
    saved=float(S["lab_momentum_edges_GeV"][-1])
    loss=-info["net_upper_GeV"]
    C=M*M-m*m-mp*mp
    out=[]
    for pmax in (.50,.60,1.,2.):
        assert pmax>saved
        e=np.hypot(mp,pmax)
        # A scalar pair can supply total energy at most 2e. The other
        # outgoing fermion then has energy no greater than 2e-m.
        fmax=2*e-m
        upper_s=(2*e)**2
        ps=np.sqrt(upper_s/4-m*m)
        amp=2*y**4*((upper_s/4)*(2*ps/C)**2+(ps*(4*m+2*M)/C)**2)
        integral=(fmax-m)*np.sqrt((fmax-m)*(fmax+m))
        G=amp*integral/(128*np.pi**3*m)
        # At least one scalar is exterior. One-exterior and two-exterior
        # configurations are disjoint, so the product's upper bound is
        # max(beta*bmax,beta**2), not their sum.
        r=loss/G
        beta=float(min(r/bmax,np.sqrt(r)))
        out.append(dict(exterior_support_max_GeV=pmax,strict_sufficient_per_mode_sup_bound=beta,
                        exact_boundary_growth_bound_check=float(G*max(beta*bmax,beta*beta)/loss)))
    return dict(scope="Frozen scalar-pair sign at tagged f=1/2 for fixed retained data. Exterior fermions can be arbitrary Pauli occupations; scalar tails require the stated finite support and pointwise upper bounds.",
                snapshot_sha256=digest(S),T_GeV=S["T_visible_GeV"],retained_outer_edge_GeV=saved,
                retained_max_scalar_occupation=bmax,protected_negative_rate_margin_GeV=loss,controls=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("paths",nargs="+",type=Path)
    a=p.parse_args()
    print(json.dumps([allowance(q) for q in a.paths],indent=2,allow_nan=False))
