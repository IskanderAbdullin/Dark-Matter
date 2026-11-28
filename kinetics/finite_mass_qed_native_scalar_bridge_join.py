"""Assemble four accepted native q subranges; do not fill absent intervals."""
import argparse
import json
from pathlib import Path

import numpy as np


BANDS=((1e-6,.0003),(.0003,.0015),(.0015,.004),(.004,.015))
LABELS=('scalar_number','scalar_energy','charged_energy')


def field(n,w,m):
    # Rows are number, scalar energy, scalar excess; columns are A,B,D.
    z=np.asarray(n,dtype=np.longdouble)
    h=(w-m)*z[0]-z[2]
    direct=w*z[0]-z[1]
    if np.max(abs(direct-h))>4e-13*np.max(abs(w*z[0])+abs(z[1])):
        raise ArithmeticError('Different scalar mass, pole, or three-sector convention')
    if np.max(abs(z[:,0]-z[:,1]-z[:,2]))>5e-13*np.max(abs(z)):
        raise ArithmeticError('Affine columns did not close')
    return np.array([z[0],z[1],h])


def roots(c):
    return [float(B/D) if B*D>0 else None for _,B,D in c]


def ledger(n,w,m):
    matrix=np.sum(n,axis=0)
    got=field(matrix,w,m)
    assembled=np.sum([field(z,w,m) for z in n],axis=0)
    if np.max(abs(got-assembled))>3e-14*np.max(abs(got)):
        raise ArithmeticError('Joining the q intervals changed the source')
    tabs=[field(z,w,m) for z in n]
    cond=np.sum(abs(np.stack(tabs)[:,:,2]),axis=0)/np.maximum(abs(got[:,2]),1e-60)
    actual=roots(assembled)
    # Keep the charge ledger as a separate linear functional. At scalar
    # number neutrality, the photon-energy flow is exactly zero.
    test=[]
    for x in sorted(set(t for t in actual if t is not None)):
        fl=np.longdouble(x)*assembled[:,2]-assembled[:,1]
        stores=np.array([fl[1],-w*fl[0],fl[2]])
        if abs(np.sum(stores))>2e-12*max(abs(stores)):
            raise ArithmeticError('Three-sector conservation failed')
        test.append(dict(n=x,scalar_photon_and_charged_energy=stores.astype(float).tolist()))
    return dict(conditional_roots=actual,slope_cancellation_numbers=cond.astype(float).tolist(),
                joined_coefficients=assembled.astype(float).tolist(),
                individual_band_roots=[roots(z) for z in tabs],tests=test)


def normalize(p,m):
    z=np.asarray(p,dtype=np.longdouble)
    if z.shape==(2,3):
        # Historical shoulder card did not retain its centered third row.
        z=np.vstack([z,z[1]-m*z[0]])
    if z.shape!=(3,3) or not np.all(np.isfinite(z)) or np.any(z[:,:2]<0):
        raise ArithmeticError('Invalid joined scalar source')
    return z


def run(root):
    def load(name):
        d=json.loads((root/name).read_text())
        if len(d)!=2 or {j['grid'] for j in d}!={'fine','third'}:
            raise ArithmeticError('An input does not contain both archive grids')
        return {r['grid']:r for r in d}
    old=load('finite-mass-qed-fixed-pole-cold-to-shoulder-gap-2026-11-28.json')
    t1=load('finite-mass-qed-native-past-bridge-2026-11-28.json')
    t2=load('finite-mass-qed-native-transfer-bridge-2026-11-28.json')
    high=load('finite-mass-qed-first-future-scalar-band-2026-11-28.json')
    result=[]
    for g in ('fine','third'):
        lo=old[g]
        if lo['joined_support']!=list(BANDS[0]) or t1[g]['scope']!=list(BANDS[1]) or t2[g]['scope']!=list(BANDS[2]):
            raise ArithmeticError('The q union no longer has exact declared support')
        w=np.longdouble(t1[g]['omega'])
        m=np.longdouble(t1[g]['m'])
        if t2[g]['omega']!=w or abs(t2[g]['m']/m-1)>2e-15:
            raise ArithmeticError('Bridge files are not for the same frozen pole')
        low=np.sum([normalize(p['raw'],m) for p in lo['pieces']],axis=0,dtype=np.longdouble)
        a,b=[np.array(t[g]['number_scalar_energy_and_centered_coefficients'],dtype=np.longdouble) for t in (t1,t2)]
        upper=normalize(high[g]['coefficients'],m)
        t=high[g]['comparative_sources'][0]
        if not t['future'][0] or abs(-t['future_photon_energy']/t['future'][0]/w-1)>2e-14:
            raise ArithmeticError('High-band photon-energy ledger has a different pole')
        for z,br in [(a,t1[g]),(b,t2[g])]:
            q=field(z,w,m)[2]
            if np.max(abs(q-br['charged_energy_coefficients']))>2e-13*np.max(abs(q)):
                raise ArithmeticError('Independent charged bridge receipt disagrees')
        result.append(dict(grid=g,pole_omega=float(w),scalar_mass=float(m),q_support=BANDS,
                           conditional_root_order=LABELS,**ledger([low,a,b,upper],w,m),
                           former_disjoint_measurement_roots=ledger([low,upper],w,m)['conditional_roots']))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path,default=Path(__file__).parent/'results',nargs='?')
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
