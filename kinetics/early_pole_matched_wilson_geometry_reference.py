"""Pole-centered collinear finite-current sensitivity at one angular coherence time."""
import argparse
import json
import numpy as np

import early_pole_finite_resonance_aperture_reference as z
from early_pole_finite_wilson_htl_reference import fold
from early_pole_finite_wilson_poles_reference import pole_fold
from early_pole_moving_soft_dipole_reference import weights


def inputs(g):
    mp, mb = np.hypot(z.X,g/2), np.hypot(z.RB*z.X,g/2)
    e,l,h,_,_=z.intervals(mp,z.MS,mb,z.Q)
    ep=(l+h)/2
    p=np.sqrt((ep-mp)*(ep+mp))
    eb=ep-e
    mu=(p*p+z.Q*z.Q+mb*mb-eb*eb)/(2*p*z.Q)
    a=z.Q*p/eb
    mD=g*np.sqrt(2/3)
    vP=p/ep
    vB=(p-z.Q*mu)/eb
    vT=z.Q*np.sqrt(1-mu*mu)/eb
    return dict(g=g,t=mD/a,vP=vP,vB_parallel=vB,vB_transverse_abs=vT,
                r=(vB-vP)*mD/a,angular_halfspan=a)


def run(level):
    grids=(4,30,18,8,30) if level==0 else (6,45,28,14,52)
    rows=[]
    for g in (.001,.003,.005):
        d=inputs(g)
        t,r,v=d['t'],d['r'],d['vP']
        for k in (8.,20.,40.):
            c=fold('triangle',t,r,v,*grids[:3],k)
            p=pole_fold(t,r,v,'triangle',*grids[3:],k)
            ref=weights(v,r,zmax=k)
            vals=dict(cut_L=c[0],cut_T=c[1],pole_L=p[0],pole_T=p[1],free_T=p[2],
                      local_cut_L=ref[1],local_cut_T=ref[2])
            # Normalized by matching the static classical HTL exponent.
            chi=(g*g/(4*np.pi*d['angular_halfspan']))*(np.sum(c)+sum(p[:2]))
            rows.append(dict(kmax=k,**d,**vals,unsubtracted_soft_exponent=chi))
    return dict(scope='Straight closed spatially collinear paths, ideal equilibrium HTL, '
                      'prescribed contour; exponent is not a production or relaxation rate.',
                level=level,rows=rows)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--level',type=int,default=0)
    print(json.dumps(run(p.parse_args().level),indent=2,allow_nan=False))
