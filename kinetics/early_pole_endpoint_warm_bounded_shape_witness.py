"""Test the same bounded in-cell perturbations for two frozen stress probes."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_stress_conditional_hull import conditioned,values


def bands(mean,width,spread):
    if not 0.<mean<1. or not 0.<width<.5:
        raise ArithmeticError('A bounded witness needs an interior normalized mean')
    h=min(width,mean,1-mean)
    if spread:
        endpoints=np.array([[0.,h],[1-h,1.]])
        t=(mean-h/2)/(1-h)
        w=np.array([1-t,t])
    else:
        endpoints=np.array([[mean-h/2,mean+h/2]])
        w=np.ones(1)
    if (np.any(w<=0) or abs(w.sum()-1)>2e-15
            or abs(w@endpoints.mean(axis=1)-mean)>3e-15):
        raise ArithmeticError('Bounded witness does not have the prescribed energy mean')
    return endpoints,w


def total(ledger,m,m0,width,spread,n):
    z,w=roots_legendre(n)
    vals=[]
    for a,b,N,mean,_,_,warm in ledger:
        if N==0:
            continue
        if warm:
            E=np.hypot(a,m0)
            dh=(b-a)*(b+a)/(np.hypot(b,m0)+E)
            interval,weight=bands(mean,width,spread)
            u=interval[:,0,None]+np.diff(interval,axis=1)*(1+z)/2
            q2=a*a+(dh*u)*(2*E+dh*u)
            P=np.einsum('j,ijk->ik',w/2,values(q2,m))
            vals.append(N*(weight@P))
        else:
            q=a+(b-a)*(1+z)/2
            pw=w*q*q
            vals.append(N*((pw@values(q*q,m))/fsum(pw)))
    out=np.array([fsum(v) for v in zip(*vals,strict=True)])
    if min(out)<=0:
        raise ArithmeticError('Expected positive conditional stresses')
    return out


def run(native,charts):
    d,y=unpack(native,False,.05)
    m0=d.state(.05,y)[0]
    tags,ledger=zip(*[conditioned(d,c,m0) for c in charts],strict=True)
    if list(tags)!=[3,4]:
        raise ArithmeticError('Need the archived G3 and G4 states')
    rows=[]
    for width in (.05,.10,.20):
        result=[]
        for M in (0.,m0,.1,1.,3.,10.,100.):
            a,b=[np.asarray([[total(c,M,m0,width,s,n) for s in (False,True)] for c in ledger])
                 for n in (48,96)]
            err=float(np.max(abs(a/b-1)))
            if err>3e-12:
                raise ArithmeticError('Resolved positive witness requires both quadratures')
            # The exact same two density families are used for both probes
            # and each test mass; no kernel-specific optimization is involved.
            row=dict(formal_mass=M,G3_G4_center_and_spread_totals=b.tolist(),
                     G4_center_vs_G3_spread_minus_one=(b[1,0]/b[0,1]-1).tolist(),
                     G4_spread_vs_G3_center_minus_one=(b[1,1]/b[0,0]-1).tolist(),
                     max_paired_relative_defect=err)
            result.append(row)
        rows.append(dict(width_cap_in_normalized_energy=width,mass_rows=result))
    return dict(frozen_x=.05,reference_mass=m0,field_order=['pressure','pseudo_pressure'],cases=rows,
                scope='Bounded piecewise densities in the reference-energy coordinate, a single matched-number/energy family for both stress probes. Formal observation masses only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--gauss',nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),[read(f) for f in a.gauss]),indent=2,allow_nan=False))
