"""Finite-mediator 5D nonthermal WEAK rate with energy-guided p segments.

No full physical Γelastic or 86eV relic conclusion. Incoming and
outgoing angular grids require separate convergence controls.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from kinetics.independent_chi_lepton_5d_quadrature import evaluate,radial_nodes


def radial_volume_error(mass=.1,cutoff=12.,segments=12,nlocal=5,
                        scheme="segmented_p"):
    p,w=radial_nodes(mass,cutoff,nlocal,scheme=scheme,segments=segments)
    exact=((mass+cutoff)**2-mass**2)**1.5/3
    return float(np.dot(w,p*p)/exact-1)


def compare(shape="oscillatory",*,cutoff=12.,max_dE=1.,nlocal=5,
            nmu=8,nc=32,nphi=16):
    nseg=math.ceil(cutoff/max_dE)
    kw=dict(nr=nlocal,nmu=nmu,nc=nc,nphi=nphi,cutoff=cutoff,
            radial_scheme="segmented_p",segments=nseg)
    r=evaluate(shape,**kw)
    r["nsegments"]=nseg
    r["radial_points_per_incoming_species"]=nseg*nlocal
    r["radial_volume_error_chi"]=radial_volume_error(.1,cutoff,nseg,nlocal)
    r["radial_volume_error_lepton"]=radial_volume_error(0.,cutoff,nseg,nlocal)
    r["max_segment_energy_width"]=cutoff/nseg
    return r


def audit(full=False):
    out={"scope":"synthetic nonthermal finite-vector Pauli incoming 5D WEAK Q/S; no cosmological Gamma/H or 86eV",
         "segmented_oscillatory":[],"segmented_smooth":[]}
    for shape in ("oscillatory","smooth"):
        for cutoff in ((8.,12.,16.,20.) if full else (8.,12.,16.)):
            out["segmented_"+shape].append(compare(shape,cutoff=cutoff))
    if full:
        out["local_refinement"]=[compare("oscillatory",nlocal=n)
                                 for n in (3,4,5,6)]
        out["angular_refinement"]=[compare("oscillatory",nmu=8,nc=nc,nphi=16)
                                   for nc in (24,32,48,64)]
    out["error_certificate"]=False
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--full",action="store_true")
    args=p.parse_args()
    print(json.dumps(audit(full=args.full),indent=2))
