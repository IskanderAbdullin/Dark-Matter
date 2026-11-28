"""Audit extra even-momentum constraints that restore heavy stress signs."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_stress_conditional_hull import conditioned,values


def local(a,b,q,w,M,power,n):
    """Taylor range certifying an observable when its finite moments are held."""
    s=q*q
    mid=(a*a+b*b)/2
    den=mid+M*M
    exponent=power-.5
    rho=(b*b-a*a)/(2*den)
    if not 0<rho<1:
        raise ArithmeticError('Heavy moment expansion escaped its cell domain')
    u=(s-mid)/den
    if abs(u).max()>rho+2e-15:
        raise ArithmeticError('Some positive mass is outside the declared support')
    w=w*s**power
    c=1.
    terms=[]
    for j in range(n+1):
        terms.append(c*fsum(w*u**j))
        c*=-((exponent+j)/(j+1))
    scale=1/(3*den**exponent)
    estimate=scale*fsum(terms)
    remainder=scale*fsum(w)*abs(c)*rho**(n+1)/(1-rho)**(exponent+n+1)
    return estimate,remainder,fsum(w)/(3*(b*b+M*M)**exponent),fsum(w)/(3*(a*a+M*M)**exponent)


def sum_rule(ledger,M,n):
    z,w=roots_legendre(96)
    rows=[]; rem=[]; direct=[]; caps=[]
    for a,b,N,mu,q,pop,warm in ledger:
        if warm:
            Y=np.asarray([local(a,b,q,pop,M,j,n) for j in (1,2)])
            rows.append(Y[:,0]); rem.append(Y[:,1]); direct.append(pop@values(q*q,M)); caps.append(Y[:,2:].T)
        else:
            x=a+(b-a)*(z+1)/2
            P=w*x*x
            v=N*(P@values(x*x,M)/fsum(P))
            rows.append(v); rem.append(np.zeros(2)); direct.append(v); caps.append(np.stack((v,v)))
    yy,err,real=[np.array([fsum(r) for r in zip(*x,strict=True)]) for x in (rows,rem,direct)]
    C=np.asarray(caps).sum(axis=0)
    limits=np.stack((np.maximum(C[0],yy-err),np.minimum(C[1],yy+err)))
    if np.any(real<limits[0]-8e-15*real) or np.any(real>limits[1]+8e-15*real):
        raise ArithmeticError('Independent nodal witness escaped analytic moment bound')
    return limits,real


def run(native,gauss):
    d,y=unpack(native,False,.05)
    m0=d.state(.05,y)[0]
    A,B=[conditioned(d,g,m0) for g in gauss]
    if [A[0],B[0]]!=[3,4]:
        raise ArithmeticError('Need actual G3 then G4')
    out=[]
    for M in (1.,3.,10.,100.):
        known=[]
        for n in range(7):
            a,c=sum_rule(A[1],M,n)
            b,d=sum_rule(B[1],M,n)
            if min(a[0])==0:
                raise ArithmeticError('Need strictly positive coarser stress lower bound')
            limits=np.array([b[0]/a[1]-1,b[1]/a[0]-1]).T
            known.append(dict(degree=n,additional_even_powers_by_field=[list(range(2*j,2*(j+n)+1,2)) for j in (1,2)],
                              conditional_G4_G3_minus_one=limits.tolist(),
                              direct_nodal_minus_one=(d/c-1).tolist()))
        out.append(dict(formal_mass=M,orders=known))
    return dict(x=.05,measurement_mass=m0,field_order=['pressure','pseudo_pressure'],rows=out,
                scope='Adds the actual independent nodal even momentum moments on each warm cell and keeps cold constant. Taylor remainder bound, ordinary floats, formal observational masses only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--gauss',nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),[read(f) for f in a.gauss]),indent=2,allow_nan=False))
