"""Locate positive pair conversion and energy exchange relative to the actual cut."""

import argparse
import json
from pathlib import Path

import numpy as np

from joint_visible_heavy_adiabatic_reference import snapshot_functions
from late_radiation_background_reference import MPL,TabulatedRadiation
from three_species_decay_onsager_reference import m,mp
from warm_coupled_serial_restart_reference import sha
from warm_light_occupation_stage_audit_reference import digest,read
from warm_pair_common_fermion_chart_reference import rate
from warm_prompt_reentry_shape_audit_reference import get_mom


def run(eos,start):
    x=read(start)
    S=x["snapshot"]
    assert x["hash"]==digest(S) and x["plan"]["eos_sha256"]==sha(eos)
    p=np.asarray(S["lab_momentum_edges_GeV"])
    C=x["plan"]["normalized_cut_GeV"]
    th=[S["entropy_variable_"+j] for j in ("scalar","fermion")]
    N,_,R=get_mom(p,th)
    T=S["T_visible_GeV"]
    H=np.sqrt((TabulatedRadiation(eos,T).state(1.)[1]+R.sum())/3.)/MPL
    B,_,F,_=snapshot_functions(S)
    controls=[]
    # Projected pair moment vector has lower exterior, represented
    # momentum shells, then upper exterior. Retain both flux directions.
    for order in ((3,7,6,2),(4,10,10,4)):
        net,out,gross,defect,_=rate(p,F,B,order)
        raw=np.stack((out[1,0],out[0,0],out[0,1],out[1,1])).reshape(2,2,2,-1)
        rows=[]
        for j,mass in enumerate((mp,m)):
            grid=np.hypot(p,mass)
            at=int(np.argmin(abs(grid-C)))
            assert abs(grid[at]-C)<2e-15 and 0<at<len(grid)-1
            lo=np.arange(raw.shape[-1])<=at
            hi=~lo
            own=np.stack([raw[j,:,0],raw[j,:,1],raw[j,:,1]+mass*raw[j,:,0]],axis=1)
            # Number normalized to that species; energies to total.
            norms=np.array([N[j].sum(),R.sum(),R.sum()])*H
            measures=[]
            for name,take in (("below",lo),("above",hi)):
                rates=own[:,:,take].sum(axis=2)/norms[None,:]
                signed=rates[0]-rates[1]
                assert min(rates[0])>=0 and min(rates[1])>=0
                measures.append(dict(side=name,gain=rates[0].tolist(),loss=rates[1].tolist(),net=signed.tolist()))
            recovered=np.array([net[j,0].sum(),net[j,1].sum(),net[j,1].sum()+mass*net[j,0].sum()])/norms
            np.testing.assert_allclose(sum(np.array(r["net"]) for r in measures),recovered,rtol=2e-12,atol=8e-15)
            rows.append(dict(species=("scalar","light")[j],rates=measures,
                             signed_total=recovered.tolist()))
        controls.append(dict(order=order,directions_pair_reaction_number_and_energy=gross.tolist(),
                             discrete_defect=float(defect),species=rows))
    return dict(input_hash=x["hash"],time=x["time"],hard_GeV=x["plan"]["hard_GeV"],energy_cut_GeV=C,
                rate_columns=["number per species N per Hubble time","kinetic energy per total rho per Hubble time",
                              "full energy per total rho per Hubble time"],
                scope="Frozen reversible pair reactions on exact daughter phases; excludes feedback and time integration.",controls=controls)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("eos",type=Path);ap.add_argument("inputs",type=Path,nargs="+")
    a=ap.parse_args()
    print(json.dumps([run(a.eos,p) for p in a.inputs],indent=2,allow_nan=False))
