"""Order limits under independent unknown normalized response outside the retained band."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_panel_moment_bounds import bounds


def fractions(a,b,top):
    """Sharp thresholds when both exterior fractions are at most epsilon."""
    a,b,top=[np.asarray(q,dtype=float) for q in (a,b,top)]
    if min(a-b)<=0. or min(b)<=0. or min(top-b)<=0. or not np.all(np.isfinite([a,b,top])):
        raise ArithmeticError('Need an already strict retained order and an exterior upper probe')
    gap=a-b
    # An adversary assigns zero probe to the lower-momentum exterior
    # of fine, and the upper exterior probe to third. Their retained
    # supports independently attain endpoints of the stated means.
    eps=gap/(gap+top)
    if not np.all((1-eps)*(-gap)+eps*top<=8*np.finfo(float).eps*top):
        raise ArithmeticError('Failed direct interpolation closure')
    return eps


def run(root,accepted):
    raw=accepted.read_bytes()
    receipt=json.loads(raw)
    m=receipt['pairs'][0]['full_pole_and_snapshot']['scalar_mass']
    options=(.0003,.0015,.015)
    out=[]
    for mesh in (1,2):
        a,b=[geometry(root/f'{g}-d2-placement-sign-endpoint-continuum-cold-{mesh}.jsonl.gz') for g in ('fine','third')]
        res={}
        for k,i,j in CHANNELS:
            x,y=[bounds(e,Z[:,i,j],m) for e,Z in (a,b)]
            if np.max(abs(y-x[:,::-1]-receipt['metrics'][k]['actual']['mean_rms_speed_shift']))>2e-12:
                raise ArithmeticError('Exterior control does not bind the accepted sign experiment')
            if not np.all(y[:,1]<x[:,0]):
                raise ArithmeticError('Observed band does not fix the ordering')
            p,q=x[:,0].copy(),y[:,1].copy()
            # RMS ordering is equivalent to second moment ordering.
            p[1]*=p[1]; q[1]*=q[1]
            z=[]
            for Q in options:
                top=np.array([Q,Q*Q,Q/np.hypot(Q,m)])
                z.append(dict(outer_q_cap=Q,max_equal_cap_per_chart=fractions(p,q,top).tolist()))
            z.append(dict(outer_q_cap=None,max_equal_cap_per_chart=[0.,0.,float((p[2]-q[2])/(1+p[2]-q[2]))]))
            res[k]=z
        out.append(res)
    for k in out[0]:
        for a,b in zip(out[0][k],out[1][k],strict=True):
            if a['outer_q_cap']!=b['outer_q_cap'] or np.max(abs(np.asarray(a['max_equal_cap_per_chart'])-b['max_equal_cap_per_chart']))>2e-9:
                raise ArithmeticError('Two meshes do not reproduce the contamination threshold')
    return dict(source_sha256=hashlib.sha256(raw).hexdigest(),columns=['mean_q','RMS_q','mean_speed_at_common_mass'],mass=m,
                response_thresholds=out[1],scope='Each chart may place an independent unknown fraction of its total normalized response outside the retained support. Thresholds are conditional sharp suprema from endpoint witnesses, not estimates of how much exterior response exists.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--root',type=Path,required=True)
    a.add_argument('--accepted',type=Path,required=True)
    p=a.parse_args()
    print(json.dumps(run(p.root,p.accepted),indent=2,allow_nan=False))
