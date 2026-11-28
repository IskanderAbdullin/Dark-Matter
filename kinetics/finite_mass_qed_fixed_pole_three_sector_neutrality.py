"""Audit the three sector energy balances on two measured disjoint scalar bands."""
import json
from math import fsum
from pathlib import Path

import numpy as np


def fields(matrix,frequency):
    t=np.array(matrix,dtype=np.longdouble)
    inv=t[:,1]
    drift=t[:,2]
    if np.min(inv)<=0 or np.min(drift)<=0:
        raise ArithmeticError('These measured cases need positive scalar inverse and drift')
    D=frequency*drift[0]-drift[1]
    B=frequency*inv[0]-inv[1]
    if D==0:
        raise ArithmeticError('Charged root unresolved at zero slope')
    return np.array([B/D,inv[1]/drift[1],inv[0]/drift[0]]),D


def flow(t,w,n):
    N,E=n*t[:,2]-t[:,1]
    G=-w*N
    C=w*N-E
    return np.array([E,G,C])


def run(root):
    def load(name):
        a=json.loads((root/name).read_text())
        if len(a)!=2 or {x['grid'] for x in a}!={'fine','third'}:
            raise ArithmeticError('Mismatched nominal scalar band archives')
        return {u['grid']:u for u in a}
    cold=load('finite-mass-qed-fixed-pole-cold-scalar-band-2026-11-28.json')
    low=load('finite-mass-qed-fixed-pole-cold-to-shoulder-gap-2026-11-28.json')
    high=load('finite-mass-qed-first-future-scalar-band-2026-11-28.json')
    out=[]
    for g in ('fine','third'):
        if low[g]['joined_support']!=[1e-6,.0003]:
            raise ArithmeticError('Low measured union does not have exact endpoints')
        w=np.longdouble(cold[g]['pole_omega'])
        L=np.asarray(low[g]['raw'],dtype=np.longdouble)
        H=np.asarray(high[g]['coefficients'],dtype=np.longdouble)[:2]
        lh,ld=fields(L,w)
        hh,hd=fields(H,w)
        both,bd=fields(L+H,w)
        if not hd<0<ld or not np.all(0<hh) or not np.all(np.diff(hh)>0):
            raise ArithmeticError('Three-sector high-q ordering did not resolve')
        R=flow(H,w,hh[0])
        if abs(fsum(R))>2e-15*max(abs(R)) or abs(R[-1])>1e-14*max(abs(R)):
            raise ArithmeticError(('Charged balance residual',g,R))
        states=[]
        for n,sg in zip((hh[0]/2,(hh[0]+hh[1])/2,(hh[1]+hh[2])/2,hh[2]*1.1),
                         ((-1,1,1),(-1,1,-1),(1,1,-1),(1,-1,-1)),strict=True):
            flux=flow(H,w,n)
            if tuple(np.sign(flux))!=sg or abs(fsum(flux))>2e-15*max(abs(flux)):
                raise ArithmeticError(('Sector sign or conservation failure',g,n,flux))
            states.append(dict(trial=float(n),scalar_photon_charged_energy=flux.astype(float).tolist()))
        r=both/hh-1
        if max(abs(r))>2e-7:
            raise ArithmeticError('Unexpected large measured low-band correction')
        out.append(dict(grid=g,conditional_root_order=['charged_energy','scalar_energy','photon_energy_and_scalar_number'],
                        high_roots=hh.astype(float).tolist(),low_roots=lh.astype(float).tolist(),
                        combined_measured_roots=both.astype(float).tolist(),signed_relative_root_change=r.astype(float).tolist(),
                        charged_energy_slope_low_high=[float(ld),float(hd)],
                        energy_balance_at_high_charged_neutral=R.astype(float).tolist(),tests=states,
                        measured_regions=[[1e-6,.0003],[.004,.015]]))
    return out


if __name__=='__main__':
    print(json.dumps(run(Path(__file__).parent/'results'),indent=2,allow_nan=False))
