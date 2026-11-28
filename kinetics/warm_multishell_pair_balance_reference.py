"""Thermal gain-loss cancellation on narrow cells of a truncated shell grid."""

import argparse
import json
import numpy as np

from archived_scalar_pair_reaction_reference import gross, m, mp
from joint_visible_heavy_adiabatic_reference import snapshot_functions
from scalar_pair_expanding_continuation_reference import rate
from warm_cold_light_pair_partition_reference import split_rate


def thermal(T, mu, pmax):
    p = np.r_[0., np.geomspace(.000035, .002, 10), np.geomspace(.0025, pmax, 22)]
    th = []
    for mass in (mp,m):
        e = np.hypot(p,mass)
        th.append(np.vstack([(e[:-1]-mu)/T, np.diff(e)/T]))
    S = dict(lab_momentum_edges_GeV=p, entropy_variable_scalar=th[0], entropy_variable_fermion=th[1])
    B,_,F,_ = snapshot_functions(S)
    out = []
    for order in ((3,7,6,2),(4,10,10,4),(6,14,14,6)):
        raw,_,_ = rate(p,th,order)
        _,_,loss = gross(p,m,mp,F,B,False,*order,project=True)
        repair,_,_,row = split_rate(p,th,order,(6,12,9,12,12),cut_cells=9)
        a = np.abs(loss[0,0,1:10])
        after = repair[1,0,1:10]/a
        assert max(abs(after)) < 2e-6
        out.append(dict(order=order, uncorrected_relative_defect=(raw[1,0,1:10]/a).tolist(),
                        repaired_relative_defect=after.tolist(), **row))
    assert max(abs(x) for x in out[-1]["repaired_relative_defect"]) < 3e-10
    return dict(T=T, mu=mu, pmax=pmax, rows=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--single", type=int, choices=(0,1))
    a = p.parse_args()
    inputs = [(0.015,-0.02,.81),(.023,-.028,.9)]
    if a.single is not None:
        inputs = [inputs[a.single]]
    print(json.dumps([thermal(*v) for v in inputs],indent=2,allow_nan=False))
