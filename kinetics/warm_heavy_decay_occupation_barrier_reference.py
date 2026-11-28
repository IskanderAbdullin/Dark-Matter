"""Pointwise heavy occupation barrier for the specified decay and visible source."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import expit,kve

from joint_visible_heavy_adiabatic_reference import parent_branch,snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp,Es,width
from warm_archived_bose_pair_weighting_reference import snapshot
from warm_late_visible_recycling_envelope_reference import COUP2,VECTOR_MASS,VECTOR_GAMMA


def thermal_injection_upper(T,delta,Kpi):
    V=VECTOR_MASS
    w=V*VECTOR_GAMMA
    s=np.hypot(V*V,w)
    # All lepton flavors, < 3 incoming mass factor over the spin average.
    vmax=3*COUP2*s/(4*np.pi*((s-V*V)**2+w*w))
    MB=M*M*T/np.pi**2 * np.exp(-M/T)*kve(2,M/T)
    S=(vmax+Kpi)*MB*np.exp(-M/T)
    margin=width*delta*(1-delta)
    # E exp(-E/T) decreases above M; this ratio holds at any momentum.
    assert T<M and S<margin
    return dict(lepton_all_s_sigma_v_upper_GeVminus2=vmax,
                pion_sigma_v_upper_GeVminus2=Kpi,MB_single_charge_density_GeV3=MB,
                max_visible_source_over_decay_margin=float(S/margin),
                remaining_extra_source_budget_in_units_Gamma_M_over_E=float(
                    delta*(1-delta)-S/width))


def frozen(s,eos,grid,angle):
    phi,sp,ferm,fp=snapshot_functions(s)
    th=np.asarray(s["entropy_variable_fermion"])
    peak=float(expit(-min(np.min(th[0]),np.min(th[0]+th[1]))))
    E=np.hypot(grid,M)
    r=parent_branch(E,phi,ferm,sp,fp,angle=angle)
    restphi=float(phi(np.array([Es]))[0])
    restferm=float(ferm(np.array([M-Es]))[0])
    f0=restphi*restferm/(1+restphi-restferm)
    f=np.r_[f0,r[:,0]]
    assert np.min(f)>=0 and np.max(f)<=peak
    T=s["T_visible_GeV"]
    H=np.sqrt(TabulatedRadiation(eos,T).state(1.)[1]/3)/MPL
    return dict(angle_nodes=angle,light_retained_peak=peak,
                max_decay_only_stationary_parent_occupation_on_probed_grid=float(max(f)),
                parent_momentum_GeV_at_max=float(np.r_[0,grid][np.argmax(f)]),
                smallest_decay_response_over_nominal_visible_H=float(min(width*M*r[:,-1]/E/H)),
                stationary_at_rest=f0)


def run(eos,pion,paths):
    pi=json.loads(Path(pion).read_bytes())["cases"]
    Kpi=max(z["controls"][-1]["pion_pair_sigma_v_cap_with_stated_support"] for z in pi)
    rows=[]
    for path in paths:
        s=snapshot(path)
        p=np.geomspace(.001,5.,100)
        trials=[frozen(s,eos,p,n) for n in (9,17)]
        z,y=[r["max_decay_only_stationary_parent_occupation_on_probed_grid"] for r in trials]
        assert abs(z-y)<1e-11
        rows.append(dict(name=Path(path).name,T_GeV=s["T_visible_GeV"],grid_size=101,
                         positive_probe_p_GeV=[p[0],p[-1]],checks=trials,
                         input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest()))
    anchor=max(d["T_GeV"] for d in rows)
    return dict(scope="Unpolarized homogeneous heavy fermion mode. If light occupancy remains "
                      "below a prescribed cap, decay/inverse decay preserves that heavy cap "
                      "even for arbitrary scalar occupancy. Adds the specified visible vacuum "
                      "sources conservatively. Additional internal heavy production/scattering "
                      "requires a separate mode-resolved upper (not supplied by archived counts).",
                parent_proper_decay_width_GeV=width,anchor_T_GeV=anchor,
                visible_source_barrier_checks=[
                    dict(cap=d,**thermal_injection_upper(anchor,d,Kpi)) for d in (.001,.1,.5)],
                frozen_cases=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("pions")
    p.add_argument("snapshots",nargs="+")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.pions,a.snapshots),indent=2,allow_nan=False))
