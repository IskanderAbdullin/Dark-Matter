"""Locate the scalar momentum carrying selected frozen spectral observables."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack

LABELS=('number','kinetic_energy_above_rest','pressure_numerator','number_weighted_velocity_squared')
PROBS=(.1,.5,.9,.99)
CUTS=(1e-5,3e-4,.01,.1,1.,3.)


def parts(a,b,f,m,n):
    a,b,f=np.broadcast_arrays(a,b,f)
    t,w=roots_legendre(n)
    h=b-a
    q=a[...,None]+h[...,None]*(t+1)/2
    E=np.hypot(q,m)
    F=(np.ones_like(q),q*q/(E+m),q*q/E,(q/E)**2)
    z=np.stack([f*h/(4*np.pi**2)*np.sum(q*q*w*u,axis=-1) for u in F],axis=-1)
    count=f*h*(b*b+b*a+a*a)/(6*np.pi**2)
    if np.any(abs(z[...,0]-count)>8e-14*abs(count)):
        raise ArithmeticError('Analytic number integral did not match')
    z[...,0]=count
    if not np.all(np.isfinite(z)) or np.min(z)<0:
        raise ArithmeticError('Invalid cell integral')
    return z


def observe(d,y,m,n):
    q=d.ds.ed
    f=y[:d.ns]
    values=parts(q[:-1],q[1:],f,m,n)
    total=np.array([fsum(v) for v in values.T])
    if np.min(total)<=0 or abs(total[0]/fsum(d.ds.w*f)-1)>2e-14:
        raise ArithmeticError('Invalid or unclosed scalar measure')
    cum=np.cumsum(values/total,axis=0)
    if np.max(abs(cum[-1]-1))>2e-14:
        raise ArithmeticError('Unclosed weighted CDF')
    def right_tail(p):
        b=q[1:]
        return [fsum(u) for u in parts(np.minimum(np.maximum(q[:-1],p),b),b,f,m,n).T]/total
    out=np.zeros((len(LABELS),len(PROBS)))
    for j in range(len(LABELS)):
        for k,p in enumerate(PROBS):
            i=np.searchsorted(cum[:,j],p)
            start=cum[i-1,j] if i else 0.
            l,r=q[i:i+2]
            def objective(s):
                return (start-p)+parts(l,l+s*(r-l),f[i],m,n)[j]/total[j]
            s=brentq(objective,0.,1.,xtol=1e-14)
            out[j,k]=l+s*(r-l)
            if abs(right_tail(out[j,k])[j]-(1-p))>5e-13:
                raise ArithmeticError('Weighted quantile did not invert its CDF')
    return total,out,np.array([right_tail(p) for p in CUTS])


def run(base,refined,x):
    (a,y),(b,z)=unpack(base,False,x),unpack(refined,True,x)
    if any(not np.array_equal(getattr(a,k),getattr(b,k)) for k in ('qd','qp','wd','wp')):
        raise ArithmeticError('The comparison must keep the charge chart fixed')
    m=a.state(x,y)[0]
    mr=b.state(x,z)[0]
    if min(m,mr)<=0:
        raise ArithmeticError('Need positive scalar pole masses')
    records=[]
    for d,v,mass in ((a,y,m),(b,z,m),(b,z,mr)):
        low,high=[observe(d,v,mass,n) for n in (24,48)]
        err=max(np.max(abs(u-w)/np.maximum(abs(w),1e-105)) for u,w in zip(low,high,strict=True))
        if err>2e-11:
            raise ArithmeticError(('Unresolved moment, quantile, or tail fraction',err))
        s,q,t=high
        records.append(dict(total=s.tolist(),moment_quantile_q=q.tolist(),fraction_at_or_above_cut_q=t.tolist(),
                            quadrature_relative_defect=float(err)))
    o,p,r=[np.array(rec['total']) for rec in records]
    if abs(r[0]/p[0]-1)>2e-14:
        raise ArithmeticError('Mass-only control changed count')
    return dict(x=x,observables=LABELS,probabilities=PROBS,tail_cuts_q=CUTS,
                pole_masses=[float(m),float(mr)],cases=('warm_native_at_own_pole','warm_refined_at_native_pole',
                    'warm_refined_at_own_pole'),records=records,
                fixed_pole_refinement_relative=(p/o-1).tolist(),mass_only_relative=(r/p-1).tolist(),
                scope='Four frozen snapshot diagnostics of numerical piecewise-constant scalar spectra. The pressure numerator omits common constants and carries no CLASS or late-time prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native',type=Path)
    p.add_argument('refined',type=Path)
    p.add_argument('--x',type=float,choices=(.01,.02,.04,.05),default=.05)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),read(a.refined),a.x),indent=2,allow_nan=False))
