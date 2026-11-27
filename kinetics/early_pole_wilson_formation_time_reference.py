"""Same straight-pole paths at several finite times, no stationary line fit."""
import argparse
import json
import numpy as np
from early_pole_finite_wilson_htl_reference import fold
from early_pole_finite_wilson_poles_reference import pole_fold
from early_pole_matched_wilson_geometry_reference import inputs


def run(level):
    grids=(4,32,18,8,30) if level==0 else (6,48,28,15,56)
    d=inputs(.005)
    g,t0,r0,v=d['g'],d['t'],d['r'],d['vP']
    pref=g*g/(4*np.pi*d['angular_halfspan'])
    rows=[]
    for k in (20.,40.):
        for u in (.25,.5,1.,2.,4.,8.):
            t,r=u*t0,u*r0
            c=fold('triangle',t,r,v,*grids[:3],k)
            p=pole_fold(t,r,v,'triangle',*grids[3:],k)
            chi=u*pref*(np.sum(c)+sum(p[:2]))
            rows.append(dict(kmax=k,time_in_inverse_angular_halfspans=u,final_r=r,
                             exponents_spacelike=(u*pref*c).tolist(),
                             timelike_exponents=(u*pref*p[:2]).tolist(),
                             free_comparator=u*pref*p[2],exponent=chi,
                             coherence_amplitude=np.exp(-chi)))
    return dict(scope='Explicit path separation accumulated through time. '
                      'A finite contour amplitude is not a stationary rate or a scalar freeze-out result.',
                level=level,**d,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    print(json.dumps(run(p.parse_args().level),indent=2,allow_nan=False))
