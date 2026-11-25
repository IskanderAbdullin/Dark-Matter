"""Rayleigh quotient of isotropic, linearized NR boson collision operator.
Perturbation delta f=f(1+f)*psi, psi=y^2-a*y-b, y=p^2/(2mT).
Number and energy orthogonal. Not a global spectral gap or halo evolution.
Compare at equal n and thermal velocity; constant identical-boson sigma.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

def moments(d):
    def moment(k,weighted):
        def fun(p):
            y=p*p
            f=np.exp(-d-y)/(-np.expm1(-d-y))
            return 2*p**(2+2*k)*f*(1+f if weighted else 1)
        return quad(fun,0,10,epsabs=1e-12,epsrel=1e-10,points=[min(np.sqrt(d),9.999)])[0]
    M=np.array([moment(k,True) for k in range(5)])
    b,a=np.linalg.solve([[M[0],M[1]],[M[1],M[2]]],M[2:4])
    c=np.array([-b,-a,1.])
    norm=sum(c[i]*c[j]*M[i+j] for i in range(3) for j in range(3))
    n=moment(0,False)
    return n,norm,a,b

def quotient(z,order=80):
    d=-np.log(z);n,norm,a,b=moments(d)
    x,w=roots_legendre(order);u=(x+1)/2;w=w/2
    # Exponential maps resolve near-critical scales in total energy and fractions.
    sl=np.log1p(80/d)
    S=d*np.expm1(sl*u);sw=w*(d+S)*sl
    I=0.
    for s,ds in zip(S,sw):
        al=np.log1p(.5*s/d)
        aa=d/s*np.expm1(al*u);aw=w*(d/s+aa)*al
        # Only b>=a, with identical in/out and permutation symmetries.
        bl=np.log((d/s+.5)/(d/s+aa))
        bb=(d/s+aa[:,None])*np.expm1(bl[:,None]*u)+aa[:,None]
        bw=w*(d/s+bb)*bl[:,None]
        f1=1/np.expm1(d+s*aa[:,None]);f2=1/np.expm1(d+s*(1-aa[:,None]))
        f3=1/np.expm1(d+s*bb);f4=1/np.expm1(d+s*(1-bb))
        # Delta psi eliminates constants and linear term; energy conserved.
        diff=2*s*s*((bb*bb-bb)-(aa[:,None]**2-aa[:,None]))
        angle=2*np.sqrt(s*aa[:,None])
        I+=float(ds*s*s*8*np.sum(aw[:,None]*bw*angle*f1*f2*(1+f3)*(1+f4)*diff**2))
    rate=I/(8*n*norm)
    return dict(z=z,order=order,n_integral=n,psi_norm=norm,psi_a=a,psi_b=b,
        bracket_integral=I,rayleigh_over_n_sigma_vT=rate,
        orthogonality='number and energy, f(1+f) inner product')

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    z=[1e-5,.01,.1,.5,.9,.99,.9999,.99999999]
    rows=[quotient(v,k) for v in z for k in [48,96,144]]
    ref=rows[2]['rayleigh_over_n_sigma_vT']
    # density g/cm^3, v_T cm/s: classical rates. This mode is not the slowest mode.
    for r in rows:
        r['enhancement_over_dilute']=r['rayleigh_over_n_sigma_vT']/ref
        r['sigma_for_one_e_fold_10Gyr_rho1e-24_vT100kms']=1/(r['rayleigh_over_n_sigma_vT']*1e-24*1e7*10e9*365.25*86400)
        print(r)
    Path(args.output).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
