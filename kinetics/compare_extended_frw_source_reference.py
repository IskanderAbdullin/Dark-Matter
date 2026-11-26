"""Compare complete FRW portal ledgers with a prescribed-background integral."""
import argparse
import json
from pathlib import Path

import numpy as np

from sm_dark_entropy_exchange_reference import StandardModel


def compare(eos,cumulative,runs):
    s0=eos.state(cumulative["reheating_T_GeV"])[2]
    out=[]
    for label,rows in runs:
        for row in rows:
            visible=row["visible_bath_feedback"]
            Tf=visible["visible_baseline_T_same_a_GeV"]
            ref=min(cumulative["rows"],key=lambda q:abs(q["T_visible_GeV"]-Tf))
            assert abs(Tf-ref["T_visible_GeV"])<3e-10
            provided=np.array(row["first_pulse_count_and_total_energy_by_channel"])+np.array(row["comoving_source_channel_count_and_total_energy"])
            y=provided[:2,0]/s0
            e=np.array([provided[0,1],provided[1:,1].sum()])
            yd=y/ref["yield_per_channel"]-1.
            ed=e/ref["net_injected_a4_energy_per_pair_channel"]-1.
            assert abs(provided[2,0]-provided[1,0])<3e-13*provided[2,0]
            if row["nsteps"]>=24:
                assert max(abs(yd))<.003 and max(abs(ed))<.003
            out.append(dict(label=label,steps=row["nsteps"],nominal_final_T_GeV=Tf,actual_T_GeV=row["final_visible_temperature_GeV"],
                            known_source_count_deviation_from_prescribed_quadrature=yd.tolist(),
                            known_source_energy_deviation_from_prescribed_quadrature=ed.tolist(),
                            final_scalar_energy_fraction=row["final_scalar_over_dark_energy"],
                            retained_final_yields=row["final_counts_div_initial_visible_entropy"]))
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("cumulative")
    p.add_argument("runs",nargs="+")
    a=p.parse_args()
    c=json.loads(Path(a.cumulative).read_text())
    runs=[(Path(f).name,json.loads(Path(f).read_text())["rows"]) for f in a.runs]
    print(json.dumps(dict(rows=compare(StandardModel(a.eos),c,runs)),indent=2,allow_nan=False))
