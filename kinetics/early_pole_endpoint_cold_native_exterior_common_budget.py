"""Separate the canceled common exterior from genuinely unmatched response."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_panel_moment_bounds import bounds


def unmatched(fine_min,third_max,upper,shared):
    a,b,H=map(lambda v:np.asarray(v,dtype=float),(fine_min,third_max,upper))
    if not 0.<=shared<1. or np.min(a-b)<=0. or np.min(b)<=0. or np.min(H-a)<=0.:
        raise ArithmeticError('Require an already ordered common band and an outside probe cap')
    # Complete sharp region for the two nonnegative unmatched masses:
    #   a*eta_fine+(H-b)*eta_third < (1-shared)*(a-b).
    r=(1-shared)*(a-b)
    return np.stack((r/a,r/(H-b),r/(a+H-b)),axis=1)


def run(root,sign):
    raw=sign.read_bytes()
    fixture=json.loads(raw)
    m=fixture['pairs'][0]['full_pole_and_snapshot']['scalar_mass']
    hs=(.0003,.0015,.015)
    out=[]
    for mesh in (1,2):
        (e,A),(f,B)=[geometry(root/f'{g}-d2-placement-sign-endpoint-continuum-cold-{mesh}.jsonl.gz')
                   for g in ('fine','third')]
        all_={}
        for label,i,j in CHANNELS:
            u,v=bounds(e,A[:,i,j],m),bounds(f,B[:,i,j],m)
            want=fixture['metrics'][label]['actual']['mean_rms_speed_shift']
            if np.max(abs(v-u[:,::-1]-want))>2e-12:
                raise ArithmeticError('Measured retained native shapes do not match the exact sign receipt')
            a,b=u[:,0].copy(),v[:,1].copy()
            a[1]*=a[1]; b[1]*=b[1]
            r=[]
            for h in hs:
                F=np.array([h,h*h,h/np.hypot(h,m)])
                T=unmatched(a,b,F,0.)
                q=np.concatenate((T.ravel(),unmatched(a,b,F,.9).ravel()))
                if np.any(q<=0.) or np.max(abs(unmatched(a,b,F,.9)/T-.1))>2e-15:
                    raise ArithmeticError('The sharp region did not scale with the remaining response')
                r.append(dict(q_cap=h,retained_fine_min=a.tolist(),retained_third_max=b.tolist(),
                              when_shared_fraction_is_zero=T.tolist(),
                              when_shared_fraction_is_point_nine=unmatched(a,b,F,.9).tolist()))
            all_[label]=r
        out.append(all_)
    for k in out[0]:
        for a,b in zip(out[0][k],out[1][k],strict=True):
            if np.max(abs(np.asarray(a['when_shared_fraction_is_zero'])-b['when_shared_fraction_is_zero']))>2e-9:
                raise ArithmeticError('Both meshes must reproduce the sensitivity thresholds')
    return dict(before_source_sha256=hashlib.sha256(raw).hexdigest(),mass=m,
                probes=['mean_q','RMS_q (second moment comparison)','mean speed'],columns=[
                    'only fine unmatched','only third unmatched','same unmatched mass in both'],
                thresholds=out[1],
                scope='Exactly shared exterior mass and distribution cancels, and retained contributions are separately normalized. The boundary controls independent residual exterior measures only. No assumed or measured exterior mass.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--sign',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.sign),indent=2,allow_nan=False))
