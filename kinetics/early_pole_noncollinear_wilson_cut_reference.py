"""Exact transverse pole velocity component in the finite spacelike contour."""
import argparse
import json
import numpy as np
from early_pole_matched_wilson_geometry_reference import inputs
from early_pole_finite_wilson_htl_reference import fold,gauss


def run(naz):
    out=[]
    d=inputs(.005)
    for cut in (8.,20.):
        t,r,v=d['t'],d['r'],d['vP']
        args=('triangle',t,r,v,4,31,19,cut)
        ref=fold(*args)
        rx=t*d['vB_transverse_abs']
        for tilt in (1.,):
            a,w=gauss(naz)
            fs=[fold(*args,rx=tilt*rx,az=(x+1)*np.pi/2) for x in a]
            val=(w/2)@np.array(fs)
            out.append(dict(kmax=cut,naz=naz,rx=tilt*rx,cut_collinear=ref.tolist(),
                            cut_noncollinear=val.tolist(),
                            ratio=(val/ref).tolist()))
    return dict(scope="Finite spacelike soft test; exact starting transverse velocity; no timelike terms here.",rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--nodes',type=int,default=5)
    print(json.dumps(run(p.parse_args().nodes),indent=2,allow_nan=False))
