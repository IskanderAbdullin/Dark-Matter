"""Decompose a bounded parent/strict extension; validate independent mass masks."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from born_thermal_quark_hard_portal_reference import XEDGES, kernel, local
from disjoint_hard_thermal_measured_source_reference import HardComplement, WINDOWS
from hard_vector_nlo_overlap_audit_reference import independent
from partial_hard_vector_nlo_reheat_reference import Correction, Fold, table


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(eos, old_native, new_native, parent12, mixed, parent14, old_max=12., new_max=14.):
    if (old_max,new_max) not in ((12.,14.),(14.,16.)):
        raise ValueError("Use one tracked pair of successive parent controls.")
    edges=np.r_[XEDGES,np.arange(13.,new_max+1.)]
    lo=table(old_native,old_max)
    hi=table(new_native,edges[-1])
    nx=len(lo[0])
    assert np.array_equal(lo[0],hi[0][:nx]) and np.array_equal(lo[2],hi[2][:nx])
    for invalid in (old_max,new_max-1.):
        try:
            table(new_native,invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("strict extension escaped parent domain")
    folds=[json.loads(x.read_text()) for x in (parent12,mixed,parent14)]
    assert [x["mass_over_T_interval"][-1] for x in folds]==[old_max,old_max,new_max]
    assert all(x["parent_born_mass_panels"]==edges.tolist() for x in folds[1:])
    F=Fold(eos,(8,21,16),6,upper=.50)
    B=[HardComplement(F,tmin=.25,n=n,x_edges=edges) for n in (4,8,12)]
    Ca=Correction(F,lo,.25,8,"linear",edges)
    Cb=Correction(F,hi,.25,8,"linear",edges)
    intervals=np.unique((WINDOWS.ravel()[:,None]/[old_max,new_max]).ravel())
    zs=np.unique(np.r_[.25,.3,.323,.336,.4,intervals-1e-6,intervals+1e-6])
    checks=[]
    for T in zs:
        u12,h12=B[-2].rates(T)
        u14,h14=B[-1].rates(T)
        unsplit=local(T,kernel(12,edges))[0]
        born=float(np.max(abs(u14/unsplit-1)))
        mask=float(np.max(abs((h12-h14)/u14)))
        low_repartition=float(np.max(abs((u12-u14)/u14)))
        A=Ca.rates(T)
        D=Cb.rates(T)
        add=independent(T,hi,"linear",20,lower_x=old_max)
        w=np.stack([add.sum(axis=0),add[0]])
        defect=float(np.max(abs((D-A-w)/D[0])))
        assert max(born,mask,low_repartition)<2e-9 and defect<2e-9
        checks.append(dict(T_GeV=T,relative_full_Born_repartition_defect=born,
              fractional_mask_order_refinement=mask,strict_new_linear_strip_defect_per_full_increment=defect,
              independent_new_strip_measured_fraction=(add[1:]/w[0]).tolist()))
    roots=[]
    for a,b,c in zip(*(x["controls"] for x in folds),strict=True):
        x,y,z=[p["measured_plus_disjoint_Born_root_GeV"] for p in (a,b,c)]
        t=a["Born_T_min_GeV"]
        assert t==b["Born_T_min_GeV"]==c["Born_T_min_GeV"] and z==y<x
        ana=[]
        for mode in range(2):
            for i in range(3):
                G=[p["results"][mode]["alpha_cases"][i] for p in (a,b,c)]
                assert len(set(u["alpha_s"] for u in G))==1
                g,h,j=[u["T_R_GeV"] for u in G]
                assert j<h<g
                ana.append(dict(interpolation=a["results"][mode]["interpolation"],alpha_s=G[0]["alpha_s"],
                                conditional_roots_GeV=[g,h,j],parent_only_shift_GeV=h-g,additional_strict_shift_GeV=j-h))
        roots.append(dict(lower_T_GeV=t,bare_optimistic_roots_GeV=[x,y],cases=ana))
    refined=[]
    for T in (.30,.323,.336,.4):
        a,b,c=[z.sample(T) for z in B]
        d=[float(np.max(abs(x-y)/c)) for x,y in ((a,b),(b,c))]
        assert d[-1]<2e-9
        refined.append(dict(T_GeV=T,cumulative_relative_order_differences=d,orders=[4,8,12]))
    return dict(root_decomposition=roots,geometry_checks=checks,extended_parent_refinement=refined,
           sha256=[sha(x) for x in (eos,old_native,new_native,parent12,mixed,parent14)],
           scope="Scope and numerical support checks for a finite mass-strip extension at prescribed coupling and clock. Pure added strip equivalence uses linear interpolation; PCHIP can also re-slope the old end interval. No missing-QCD bound.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("eos","old_native","new_native","parent12","mixed","parent14"):
        p.add_argument(name,type=Path)
    p.add_argument("--old-max",type=float,default=12.)
    p.add_argument("--new-max",type=float,default=14.)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.old_native,a.new_native,a.parent12,a.mixed,a.parent14,a.old_max,a.new_max),indent=2,allow_nan=False))
