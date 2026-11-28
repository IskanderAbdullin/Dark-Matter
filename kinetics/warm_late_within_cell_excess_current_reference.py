"""Perturb one occupied cell while preserving its full and cut number/energy."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M, m, mp
from warm_corrected_exterior_reactivated_tail_stress_reference import reactivate
from warm_corrected_late_prompt_cut_aligned_reference import split_at_cut
from warm_light_occupation_stage_audit_reference import digest, read
from warm_prompt_excess_energy_tail_reference import integrate
from warm_prompt_reentry_shape_audit_reference import get_mom


def invariant_tiles(fn, lo, hi, mass, species_g, cut):
    # All three tiles are strictly inside the same original represented cell.
    edges=np.linspace(lo+.04*(hi-lo),hi-.04*(hi-lo),4)
    assert lo<edges[0]<edges[-1]<hi and lo>cut
    def fold(n):
        x,w=roots_legendre(n)
        half=np.diff(edges)[:,None]/2
        E=edges[:-1,None]+half*(1+x)
        dN=species_g*E*np.sqrt((E-mass)*(E+mass))/(2*np.pi**2)*fn(E)*(half*w)
        return np.stack([dN.sum(axis=1),(dN*(E-cut)).sum(axis=1)])
    V=fold(24)
    X=V[1]/V[0]
    coefficient=np.array([X[1]-X[2],X[2]-X[0],X[0]-X[1]])/V[0]
    coefficient/=np.max(abs(coefficient))
    independent=fold(36)
    rel=(independent @ coefficient) / (abs(independent) @ abs(coefficient))
    assert np.max(abs(rel))<5e-13
    return edges,coefficient,rel


def remap(fn, old_edges, edges, coeff, amplitude):
    def F(E):
        k=np.searchsorted(edges,E,side="right")-1
        offset=np.where((k>=0)&(k<len(coeff)),coeff[np.clip(k,0,len(coeff)-1)],0.)
        return fn(E)*(1.+amplitude*offset)
    return F,np.unique(np.r_[old_edges,edges])


def run(eos, prior, archive, reference):
    x=read(archive)
    ref=read(reference)
    S=x["snapshot"]
    assert digest(S)==x["snapshot_sha256"]==ref["source_sha256"]
    row=next(a for a in ref["controls"] if a["physical_hard_shell_halfwidth_GeV"]==.0025)
    center=row["local_crossing_candidate_GeV"]
    D,*_=reactivate(S,center,.0025,cap=1.05)
    C=float(np.hypot(M,1./read(prior)["branches"][-1]["entropy_expansion_from_original_anchor"]))
    D=split_at_cut(D,C)
    p=np.asarray(D["lab_momentum_edges_GeV"])
    E=get_mom(p,[D["entropy_variable_"+v] for v in ("scalar","fermion")])[2].sum()
    H=np.sqrt((TabulatedRadiation(eos,D["T_visible_GeV"]).state(1.)[1]+E)/3.)/MPL
    F,A,G,B=snapshot_functions(D)
    original=integrate(D,C,12,True)
    explicit=integrate(D,C,12,True,evaluators=(F,A,G,B))
    assert np.array_equal(original,explicit)
    original=original/(E*H)
    hi=center+.0025
    ix=np.flatnonzero(p==hi)
    assert len(ix)==1 and ix[0]>0
    j=ix[0]
    rows=[]
    for sp,fn,edges,mass,g in (("scalar",F,A,mp,1),("light",G,B,m,4)):
        lo=max(edges[j-1],C)+1e-11
        tiles,coeff,err=invariant_tiles(fn,lo,edges[j],mass,g,C)
        rates=[]
        for amp in (-.5,0.,.5):
            changed,X=remap(fn,edges,tiles,coeff,amp)
            span=np.linspace(tiles[0],tiles[-1],101)
            assert np.all(changed(span)>=0)
            if sp=="light":
                # f is monotone inside the base cell and multiplication
                # is constant within each tile; the extreme endpoints
                # enclose every multiplier and base occupation.
                assert 1.5*float(np.max(fn(tiles[[0,-1]])))<1.
            ev=(changed,X,G,B) if sp=="scalar" else (F,A,changed,X)
            checks=[]
            for order in (8,12,16):
                out=integrate(D,C,order,True,evaluators=ev)/(E*H)
                checks.append([float(out[3]),float(out[1]+out[2])])
            assert max(abs(r[0]-checks[-1][0]) for r in checks)<2e-8*checks[-1][1]
            rates.append(dict(amplitude=amp,orders=[8,12,16],signed_and_gross=checks))
        assert rates[0]['signed_and_gross'][-1][0]>0>rates[-1]['signed_and_gross'][-1][0]
        assert abs(rates[1]['signed_and_gross'][-1][0]-original[3])<1e-15*(original[1]+original[2])
        rows.append(dict(species=sp,base_interval_GeV=[float(edges[j-1]),float(edges[j])],
                         tile_edges_GeV=tiles.tolist(),coefficients=coeff.tolist(),
                         relative_null_errors=err.tolist(),rates=rates))
    return dict(scope="Positive held within-cell phase perturbations, zero changes of count and energy "
                      "in affected species and cell, entirely above the same cut; frozen prompt only.",
                archive_sha256=x["snapshot_sha256"],base_sha256=digest(D),
                physical_cut_GeV=C, packet_center_GeV=center,controls=rows,
                unchanged_default_angular_result=original.tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","prior","archive","reference"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.prior,a.archive,a.reference),indent=2,allow_nan=False))
