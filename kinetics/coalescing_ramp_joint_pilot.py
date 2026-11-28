"""Convolve the two-root primitive and integrate the real finite-mass spinors."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from finite_mass_qed_coalescing_histogram_edge_reference import actual, cumulative


def anti(t,h):
    def atan_int(x):
        return x*np.arctan(x)-.5*np.log1p(x*x)

    def log_int(x):
        return x*np.log1p(x*x)-2*x+2*np.arctan(x)

    # At h=0 the quotient in the cumulative tends to 2t/(1+t^2).
    ell=np.log1p(t*t) if h==0 else (log_int(t)-log_int(t-h))/h
    return np.pi*t+atan_int(t)+atan_int(t-h)+ell


def analytic(w,h,eta,width):
    l,r=eta-width,eta+width
    linear=(anti(r,h)-anti(l,h))/(r-l)
    cubic=quad(lambda z:.75*(1-z*z)*cumulative(eta+z*width,h),-1.,1.,
               epsabs=2e-12,epsrel=2e-12)[0]
    return w[1]*2*np.pi+(w[0]-w[1])*np.array([linear,cubic])


def direct(w,h,eta,width,smooth):
    def integrand(x):
        z=max(0.,min(1.,.5+(x-eta)/(2*width)))
        if smooth:
            z=z*z*(3-2*z)
        f=(h*h+4)/((x*x+1)*((h-x)**2+1))
        return (w[0]+(w[1]-w[0])*z)*f

    edges=(-np.inf,eta-width,eta+width,np.inf)
    return sum(quad(integrand,a,b,epsabs=2e-12,epsrel=2e-12)[0]
               for a,b in zip(edges[:-1],edges[1:],strict=True))


def run():
    results=Path(__file__).parent/'results'
    saved=json.loads((results/'finite-mass-qed-coalescing-histogram-edge-2026-11-28.json').read_text())['data']
    resolved=[]
    actual_rows=[]
    broad=[]
    for grid in ('fine','third'):
        w=np.array(next(x for x in saved if x['grid']==grid)['normalized_step_weights'])
        for h in (0.,.5,2.,6.,16.):
            for eta in (-2.,0.,3.):
                sharp=w[1]*2*np.pi+(w[0]-w[1])*cumulative(eta,h)
                for width in (.08,.5,1.,3.,12.):
                    sm=analytic(w,h,eta,width)
                    independently=[direct(w,h,eta,width,is_cubic) for is_cubic in (False,True)]
                    err=max(abs(sm-independently))
                    reflected=analytic(w,h,h-eta,width)
                    reflection=max(abs(sm+reflected-2*np.pi*sum(w)))
                    if max(err,reflection)>2e-10:
                        raise ArithmeticError((grid,h,eta,width,err,reflection))
                    resolved.append(dict(grid=grid,separation=h,edge_position=eta,halfwidth=width,
                                         linear=sm[0],cubic=sm[1],sharp=sharp,
                                         direct_defect=err,reflection_defect=reflection))

        for h,eta in ((0.,3.),(2.,0.),(6.,-2.),(16.,0.)):
            a=h/2-eta
            center=np.pi*sum(w)
            coefficient=np.pi*(w[1]-w[0])*a*np.array([1.,1.5])
            for sigma in (32.,128.,512.):
                sm=analytic(w,h,eta,sigma)
                relative=max(abs((sigma*(sm-center)-coefficient)/coefficient))
                if sigma==512. and relative>4e-4:
                    raise ArithmeticError(('Wide ramp first correction',grid,h,eta,relative))
                broad.append(dict(grid=grid,separation=h,edge_position=eta,halfwidth=sigma,
                                  midpoint_area=center,area=sm.tolist(),
                                  scaled_limit_correction=coefficient.tolist(),
                                  first_correction_relative_error=relative))

        with gzip.open(results/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            endpoint=json.load(f)
        widths=(.5,3.,12.)
        for h,eta in ((0.,0.),(0.,3.),(.5,-2.),(2.,0.),(6.,3.),(16.,0.)):
            for eps in (2e-10,4e-11):
                row=actual(endpoint,grid,h,eta,eps,halfwidths=widths)
                w=row['normalized_step_weights']
                target=np.array([analytic(w,h,eta,d) for d in widths])[:,:,None]
                target=target*np.array([1.,h*h/(h*h+4)])[None,None,:]
                defect=float(np.max(abs(row['smoothed_spinor_integrals']-target)))
                if defect>8e-5:
                    raise ArithmeticError(('finite-mass smoothed check failed',row,target,defect))
                row.update(smoothed_halfwidths=widths,smoothed_joint_limits=target.tolist(),
                           maximum_smoothed_finite_mass_defect=defect)
                actual_rows.append(row)

    return dict(scope='Prescribed two-edge smoothing profiles; diagnostic denominators, not a kinetic rate.',
                leading_kernel=resolved,broad_ramps=broad,finite_mass=actual_rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
