"""Quantum Bose 2<->2 detailed-balance witness for Taylor-lowered models.
Requires numpy and reproduce_cutoff_ladder.py, does NOT compute full C22.
"""
import json
import numpy as np
from reproduce_cutoff_ladder import occupation

def main(eta=4.):
    t=np.array([.05,.35,.2,.2])
    witnesses=[]
    for k in range(7):
        f=1/np.expm1(eta-t) if k==0 else occupation(t,eta,k)
        L=np.log1p(1/f)
        gap=float(L[0]+L[1]-L[2]-L[3])
        witnesses.append({'k':k,'imbalance_Delta_L':gap})
    assert abs(witnesses[0]['imbalance_Delta_L'])<1e-12
    assert all(x['imbalance_Delta_L']>.5 for x in witnesses[1:])
    return {'eta':eta,'collision_depths_t':t.tolist(),
       'energy_sum_equal':bool(abs(t[0]+t[1]-t[2]-t[3])<1e-12),
       'witnesses':witnesses,'checks_passed':2,
       'limitation':'Detailed-balance test for closed ergodic quantum elastic collisions only; no full C22 integral, rate, or open-halo escape flux.'}
if __name__=='__main__': print(json.dumps(main(),indent=2))
