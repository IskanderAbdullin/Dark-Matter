"""Equal-number conversion cannot change the scalar-equivalent count inventory."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import m, mp

R=mp/m
TARGET=.120


def read(p):
    return json.loads(p.read_text())


def entries(label, omega):
    f,c=np.asarray(omega)
    assert f>=0. and c>=0.
    y=f+R*c
    # If all stable counts became scalar rest mass this is the most
    # scalar density possible, irrespective of redistribution/kinetics.
    assert y>0.
    hi=y/R
    out=dict(case=label,snapshot_scalar_stable_Omega_h2=[f,c],scalar_equivalent_inventory=y,
             all_stable_endpoint=hi,absolute_scalar_fraction_of_target_upper=y/TARGET)
    if y<=TARGET<=hi:
        phi=(y-R*TARGET)/(1.-R)
        st=TARGET-phi
        assert 0.<=phi<=y and 0.<=st<=hi
        # Independent linear solve with a dimensionless, well-conditioned matrix.
        p=np.linalg.solve(np.array([[1.,1.],[1.,R]]),np.array([TARGET,y]))
        assert np.max(abs(p-[phi,st]))<1e-15
        out.update(closed_partition_at_target_scalar_stable=[phi,st],
                   resulting_scalar_fraction_if_total_matches_target=phi/TARGET)
    else:
        out.update(closed_partition_at_target_scalar_stable=None,resulting_scalar_fraction_if_total_matches_target=None)
    out["required_inventory_multipliers_at_target"]={
        str(q):TARGET*(q+R*(1.-q))/y for q in (.5,.9,1.)}
    return out


def run(dist,early,warm,late):
    A,B,C,D=map(read,(dist,early,warm,late))
    E=A["eos_sha256"]
    assert E==B["eos_sha256"]
    rows=[]
    residual=[]
    for i,(x,s) in enumerate(zip(A["cases"],B["episodes"],strict=True)):
        assert x["T_R_GeV"]==s["T_R_GeV"] and x["Tmax_GeV"]==s["Tmax_GeV"]
        a=np.array([b["Omega_h2"] for b in x["species"]])
        delta=float(np.asarray(s["modes"][-1]["new_scalar_omega"])[:,0].sum())
        assert 0.<delta<R*a[1]
        b=a+delta*np.array([1.,-1./R])
        old=entries(f"low_reheat_prompt_{i}",a)
        new=entries(f"low_reheat_early_forward_proxy_{i}",b)
        err=abs(old["scalar_equivalent_inventory"]/new["scalar_equivalent_inventory"]-1.)
        assert err<2e-14
        residual.append(float(err))
        rows.extend([old,new])
    aa=[]
    for z,x in zip(("warm_q16_12MeV","warm_q32_12MeV"),C["cases"],strict=True):
        aa.append(entries(z,x["inclusive_fixed_Omega_h2"]))
    point=D["chained_endpoint_cases"][-1]
    after=entries("warm_q16_1.206MeV",point["retained_and_all_particle_tail_Omega_h2"])
    drift=abs(after["scalar_equivalent_inventory"]/aa[0]["scalar_equivalent_inventory"]-1.)
    assert drift<5e-9
    rows.extend(aa+[after])
    return dict(target_total_Omega_h2=TARGET,mp_over_stable_mass=R,cases=rows,
                low_early_forward_inventory_identity_relative_defects=residual,
                warm_closed_count_12_to_1p206_relative_defect=float(drift),
                sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (dist,early,warm,late)},
                scope="Rest-count inventory comparison under fixed visible comoving entropy and no further net number injection/removal in the two tracked species. Equal-number conversion or elastic scattering alone cannot raise this inventory. This is not a reaction-rate, thermalization, transfer, observational, or existence claim for an alternative physical model.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("dist","early","warm","late"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.dist,a.early,a.warm,a.late),indent=2,allow_nan=False))
