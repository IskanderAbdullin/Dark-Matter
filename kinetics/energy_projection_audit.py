"""Test conservation of parent energy by the occupation-weighted 1<->2 operator.

A barycentric projection conserves energy for occupation-INDEPENDENT event
weights. Actual rates are evaluated separately in four recipient-bin pairs,
so energy conservation must be tested after multiplying by their own S_ij.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from one_to_two_fixed_bath import FlatThermalDecay


def diagnose(nq=100, np1=140, nmu=72, populated=True):
    k=FlatThermalDecay(nq=nq,np1=np1,nmu=nmu)
    y=np.zeros(2*nq)
    if populated:
        y[:nq]=100.0*np.exp(-k.avg_q/0.07)
        y[nq:]=0.25*np.exp(-k.avg_q)
    energy_error=0.0
    absolute_reference=0.0
    per_event=np.zeros_like(k.parentE)
    for ip,ic,rate in k.event_gain_loss(y,stimulated=True):
        diff=k.avg_q[ip]+k.avg_q[ic]-k.parentE
        energy_error+=float(np.dot(rate,diff))
        absolute_reference+=float(np.dot(np.abs(rate),k.parentE))
        per_event+=rate*diff
    return {
        "grid":[nq,np1,nmu],
        "populated":bool(populated),
        "profile":"f_phi=100 exp(-q/0.07), f_chi2=0.25 exp(-q)" if populated else "vacuum",
        "energy_discrepancy_over_T4_per_Gamma":energy_error,
        "signed_relative_energy_discrepancy":energy_error/absolute_reference,
        "absolute_per_parent_event_relative_error":float(np.abs(per_event).sum()/absolute_reference),
        "absolute_reference":absolute_reference,
        "warning":"Fixed parent reservoir, so this is a diagnostic difference between daughter energy deposition and parent energy lost in the same weighted reactions, not cosmological energy conservation."
    }


def run(quick=False):
    configs=[(50,64,48)] if quick else [(50,64,48),(100,140,72),(150,180,96),(180,240,120)]
    rows=[{"vacuum":diagnose(*c,populated=False),"populated":diagnose(*c,populated=True)} for c in configs]
    for pair in rows:
        assert abs(pair["vacuum"]["signed_relative_energy_discrepancy"])<1e-11
        assert pair["populated"]["absolute_per_parent_event_relative_error"]>0
    return {"baseline_commit":"ea6ec461997600d5f00bf244ea808607d40ef8dd","scenarios":rows}


if __name__=="__main__":
    pa=argparse.ArgumentParser()
    pa.add_argument("--quick",action="store_true")
    pa.add_argument("--output",default="outputs/kinetic_energy_projection_audit.json")
    a=pa.parse_args()
    out=run(a.quick)
    p=Path(a.output)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(out,indent=2),encoding="utf-8")
    print(json.dumps(out,indent=2))
