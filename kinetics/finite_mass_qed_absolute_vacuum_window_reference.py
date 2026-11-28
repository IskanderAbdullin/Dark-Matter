"""Vacuum MS-bar Yukawa normalization and infrared-safe rest-cut control."""
import gzip
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from finite_mass_qed_vacuum_scalar_vertex_shape_reference import vertex,radiation,vacuum_recoil

def Zfin(m,mu,photon):
    # Finite Z_pole in Feynman gauge, after removal of its common UV pole.
    l=photon**2; A=m*m
    D=lambda x:x*x*A+(1-x)*l
    w=np.sqrt(4*A*l-l*l)
    j=np.log(A/l)/(2*A)+l/(A*w)*(np.arctan2(2*A-l,w)-np.arctan2(-l,w))
    x=np.array([.2,1,4,16,64])*photon/m
    x=x[x<1]
    return (1+quad(lambda u:2*(1-u)*np.log(D(u)/mu**2),0,1,points=x,epsabs=1e-12)[0]
            +4*A*(j-quad(lambda u:u**3/D(u),0,1,points=x,epsabs=1e-12)[0]))

def asym(M,m,q2,R):
    f,J=vertex(M,m,q2)
    B,C,_=radiation(M,m,q2)
    V=f-8+3*np.log(m/M)
    recoil=vacuum_recoil(M,m,q2,R)
    ans=.5*V+B*np.log(R)+C+recoil
    rows=[]
    for l in (1e-4,1e-6,1e-8):
        ph=M*l
        vr=vertex(M,m,q2,ph)[0]-4-J*np.log(1/l)
        vr+=.5*(Zfin(M,M,ph)+Zfin(m,M,ph))
        phi=np.sqrt((R-l)*(R+l))/R
        ang=radiation(M,m,q2)[2]
        r=B*np.log(R/l)+quad(ang,0,phi,epsabs=5e-12)[0]+recoil
        rows.append(dict(photon_mass_over_M=l,combined=float(.5*vr+r),defect=float(abs(.5*vr+r-ans))))
    return ans,rows


def self_test(m,mu):
    lam=.001*m
    q=Zfin(m,mu,lam)
    with mp.workdps(66):
        m,mu,lam=map(mp.mpf,(m,mu,lam))
        def t(a,x):
            return x*m*m+(1-x)*lam*lam-x*(1-x)*a
        cuts=map(mp.mpf,('0','.0001','.001','.01','.1','1'))
        cuts=list(cuts)
        def f(a):
            return 1+mp.quad(lambda x:2*(1-x)*mp.log(t(a,x)/mu**2),cuts)
        def summ(a):
            return 1+mp.quad(lambda x:-2*(1+x)*mp.log(t(a,x)/mu**2),cuts)
        e=mp.mpf('1e-9')*lam**2
        c=(summ(m*m+e)-summ(m*m-e))/(2*e)
        check=float(abs(mp.mpf(q)-(f(m*m)+2*m*m*c)))
        if check>2e-11:
            raise ValueError('Wavefunction residue failed the self-energy derivative check.')
    return check


def run():
    M0=.05
    f,J=vertex(M0,M0,0.)
    if max(abs(f-8+2),abs(J-4))>2e-13:
        raise ValueError('The zero-transfer equal-mass scalar mass derivative failed.')
    ratio=.003/M0
    e0=self_test(M0,M0)
    e1=self_test(M0*ratio,M0)
    rows=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as fi:d=json.load(fi)
        M=d['x_end']; m=M*d['mass_ratio']; s=d['checkpoints'][-1]['scalar_pole_over_T']
        lam=M*1e-5
        # Pure UV scale response at a fixed common infrared subtraction point.
        f,J=vertex(M,m,s*s,lam)
        def scalechange(k):
            scale=M*np.exp(k)
            return 8*k+(Zfin(M,scale,lam)-Zfin(M,M,lam)+Zfin(m,scale,lam)-Zfin(m,M,lam))/2
        rescale=max(abs(scalechange(k)-6*k) for k in (-np.log(2.),np.log(2.)))
        if rescale>2e-12:
            raise ValueError('The finite scalar-current scale dependence did not reproduce.')
        rr=[]
        for r in (.001,.01,.05):
            val,ladder=asym(M,m,s*s,r)
            delta=val-asym(M,m,0.,r)[0]
            if ladder[-1]['defect']>5e-7:
                raise ValueError('The regulator-free finite window did not reproduce.')
            rr.append(dict(rest_photon_cap_over_parent_mass=r,
                           absolute_vacuum_coefficient_before_gauge=float(val),
                           coefficient_relative_to_zero_transfer=float(delta),
                           at_g_005_relative_to_born=float(.005**2/(4*np.pi**2)*val),
                           regulator_ladder=ladder))
        rows.append(dict(grid=grid,mu_over_M=1.,scheme_scale_defect=float(rescale),rows=rr))
    return dict(scope='Vacuum O(g²) on-shell fermion, MS-bar Yukawa at parent scale, parent-rest photon cap; no finite thermal contributions.',
                equal_mass_zero_transfer_scalar_factor=float(vertex(M0,M0,0.)[0]-8),
                independent_self_energy_derivative_defects=[e0,e1],rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
