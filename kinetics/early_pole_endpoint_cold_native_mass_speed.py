"""Compare cold-band mean speeds at common mass and at the two native pole masses."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_native_placement_refinement import CHANNELS, probability
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_cold_panel_moment_bounds import bounds
from early_pole_endpoint_continuum_cold_band_receipt import run as receipt


def velocities(q,m):
    return q/np.hypot(m,q)


def evaluate(a,A,b,B,m,n):
    if not (0.<n<m and max(a[-1],b[-1])<n/np.sqrt(2) and a[0]==b[0] and a[-1]==b[-1]):
        raise ArithmeticError('This monotone mass correction is restricted to a matched cold band')
    z={}
    for name,i,j in CHANNELS:
        f,g=probability(A[:,i,j]),probability(B[:,i,j])
        u=bounds(a,f,m)[-1]
        shared=bounds(b,g,m)[-1]-u[::-1]
        native=bounds(b,g,n)[-1]-u[::-1]
        shift=velocities(b,n)-velocities(b,m)
        c=[fsum(g*shift[:-1]),fsum(g*shift[1:])]
        if (np.any(np.diff(shift)<=0.) or native[0]-shared[0]<c[0]-1e-14
                or native[1]-shared[1]>c[1]+1e-14 or native[0]<shared[0] or native[1]<shared[1]):
            raise ArithmeticError('Mass correction did not match endpoint-moment bounds')
        def sign(d):
            return 'third_lower' if d[1]<-2e-12 else 'third_higher' if d[0]>2e-12 else 'unresolved'
        z[name]=dict(shared_mass_speed_difference=shared.tolist(),native_mass_speed_difference=native.tolist(),
                     third_mass_correction_only=c,signed_order=[sign(shared),sign(native)])
    return z


def run(root,archive,lock):
    data=lock.read_bytes()
    d=json.loads(data)
    poles=[x['full_pole_and_snapshot'] for x in d['pairs']]
    if [(x['grid'],x['time_x']) for x in poles]!=[('fine',.05),('third',.05)]:
        raise ArithmeticError('Expected the held native D2 endpoint comparison')
    paths=[[root/f'{g}-d2-placement-sign-endpoint-continuum-cold-{j}.jsonl.gz' for j in (1,2)] for g in ('fine','third')]
    if any(receipt(p,archive)!=prior for p,prior in zip(paths,d['pairs'],strict=True)):
        raise ArithmeticError('Exact optics and original refined cold-band ledgers did not reproduce')
    m,n=[p['scalar_mass'] for p in poles]
    rows=[]
    for k,l in zip(*paths,strict=True):
        a,A=geometry(k)
        b,B=geometry(l)
        rows.append(evaluate(a,A,b,B,m,n))
    for name,_,_ in CHANNELS:
        for key,a in rows[0][name].items():
            b=rows[1][name][key]
            if (a!=b if key=='signed_order' else np.max(abs(np.asarray(a)-b))>2e-10):
                raise ArithmeticError('Native mass attribution did not reproduce on both meshes')
    return dict(accepted_band_receipt_sha256=hashlib.sha256(data).hexdigest(),q_band=[1e-6,1.5e-5],
                fine_and_third_native_scalar_masses=[m,n],native_minus_fine_mass_relative=n/m-1,
                band_speed_results=rows[1],
                scope='Mean speed proxy for two held chart masses only. No speed ordering beyond the measured q band, kinetic history, or observable abundance inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','archive','lock'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.archive,a.lock),indent=2,allow_nan=False))
