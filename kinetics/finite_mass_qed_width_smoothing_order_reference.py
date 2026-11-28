"""Local comparison of fixed vs regulator-tracking occupation smoothing."""
from pathlib import Path
import json
import numpy as np
from scipy.special import erf,roots_legendre

from finite_mass_qed_correlated_overlap_guards_reference import polars


def model(a,b,sign,n,m,lim):
    def f(x):
        c=abs(a*x)
        y,w=polars(-1e13,1e13,b*x,max(1e-40,min(c,1.)),1.,n)
        X=np.sqrt(y*y+c*c)
        Z=np.sqrt(y*y+c*c-1j*sign)
        Δ=(-1j*sign)/(2*X*Z*(X+Z))
        return np.sum(w*Δ/(y-b*x+1j))
    br=np.geomspace(1e-10,lim,105)
    x,w=roots_legendre(m)
    left,right=np.log(br[:-1]),np.log(br[1:])
    hh=(right-left)/2
    xs=np.exp(left[:,None]+hh[:,None]*(x+1))
    vs=np.array([[f(q).real for q in a] for a in xs])*xs*hh[:,None]*w
    norm=np.sum(vs)
    if not abs(norm)>1e-9:
        raise ArithmeticError('Chosen width path has no step real coefficient.')
    rs=np.array([1e-6,.01,.1,1,10,100,1e3,1e4,1e6,1e8,1e10])
    def r(s):
        gf=erf(xs/s)
        z=np.minimum(1.,xs/s)
        g2=z*z*(3-2*z)
        return (np.sum(vs*gf)/norm,np.sum(vs*g2)/norm)
    return np.array([r(s) for s in rs]),rs


def run():
    root=Path(__file__).parent/'results'
    with open(root/'finite-mass-qed-correlated-overlap-guards-2026-11-27.json') as f:
        raw=json.load(f)
    result=[]
    for z in raw['rows']:
        row=next(x for x in z['rows'] if x['tuned_step'])
        a,b=row['local_relative_angular_scales']
        sign=-1 if row['cut']=='daughter' else 1
        high,r=model(a,b,sign,26,12,1e10)
        low,_=model(a,b,sign,19,9,1e9)
        if max(abs(high.ravel()-low.ravel()))>8e-5:
            raise ArithmeticError('Joint width and occupation smoothing did not resolve.')
        if max(abs(high[-1]))>2e-4:
            raise ArithmeticError('Fixed physical smoothing limit did not vanish.')
        result.append(dict(grid=z['grid'],cut=row['cut'],relative_smoothing_scale=r.tolist(),
                           normalized_real_regulator_offset=high.tolist(),
                           order_and_grid_defect=float(max(abs(high.ravel()-low.ravel())))))
    return dict(scope='Manufactured local smooth transitions on two charged overlap fixtures; not a spectral width or conservation remap.',
                model_order=['erf','cubic_smoothstep'],rows=result)


if __name__=='__main__': print(json.dumps(run(),indent=2,allow_nan=False))
