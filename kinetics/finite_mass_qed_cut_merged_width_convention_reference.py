"""Manufactured constant-residue convention test anchored to an actual merged ratio."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad


def run():
    root=Path(__file__).parent/'results'
    states=json.loads((root/'finite-mass-qed-cut-pole-coalescence-2026-11-27.json').read_text())['data']
    ratios=[v['coherent_residue_over_incoherent_diagonals'] for s in states for x in s['checks'] for v in x['rows']]
    rows=[]
    for r in (min(ratios),max(ratios)):
        # Two real orthogonal spin components reproduce the exact total
        # diagonal norm 1 and the measured difference-squared ratio r.
        a=np.array([np.sqrt(1.-r/2),np.sqrt(r/2)])/np.sqrt(2)
        b=np.array([a[0],-a[1]])
        dim=[]
        for h in (0.,.01,.1,1.,10.,60.):
            spec=[]
            for sign in (1,-1):
                def f(u):
                    amp=a/(u+1j)+b/(h-u+1j*sign)
                    return abs(amp)@abs(amp)/np.pi
                p0,p1,p2=(quad(f,l,u,epsabs=1e-17,epsrel=2e-10,limit=360)[0]
                          for l,u in ((-np.inf,0.),(0.,h),(h,np.inf)))
                exact=1. if sign==1 else (r+h*h/4)/(1.+h*h/4)
                absolute=abs(p0+p1+p2-exact)
                if absolute>1e-9*exact+1e-17:
                    raise ArithmeticError(('Constant-residue width test failed',r,h,sign,absolute))
                spec.append(dict(second_imaginary_sign=sign,scaled_diagonal_subtracted_analytic=exact-1.,
                                scaled_unsubtracted_analytic=exact,independent_abs_defect=absolute))
            dim.append(dict(invariant_separation_over_first_width=h,conventions=spec))
        rows.append(dict(anchored_coherent_residue_ratio=r,dimensionless_width_fixture=dim))
    return dict(scope='Constant-residue algebra, comparing two manufactured regulator sign choices only; neither chooses the physical contour or finite matching.',data=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
