"""Unequal-mass vacuum scalar vertex shape after a common soft subtraction."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from early_pole_hard_yukawa_photoproduction_reference import slash, MET
from finite_mass_qed_signed_real_soft_recoil_reference import setup, angular_scaled


def contraction_check(M,m,q2):
    P=np.array([M,0.,0.,0.])
    Eb=(M*M+m*m-q2)/(2*M)
    b=np.array([Eb,0.,0.,np.sqrt((Eb-m)*(Eb+m))])
    G=slash(np.eye(4))
    I=np.eye(4,dtype=complex)
    U=slash(b)+m*I
    V=slash(P)+M*I
    denominator=np.trace(U@V)
    d=M*Eb
    out=[]
    for u,t in ((.03,.11),(.18,.45),(.37,.8),(.6,.19),(.91,.9)):
        x,y=u*t,u*(1-t)
        L=x*P+y*b
        A=slash(b-L)+m*I
        B=slash(P-L)+M*I
        mat=sum(g*Z@A@B@Z for g,Z in zip(MET,G,strict=True))
        actual=np.trace(U@mat@V)/denominator
        D=t*t*M*M+(1-t)**2*m*m+2*t*(1-t)*d
        tr=t*(M*M+d)+(1-t)*(m*m+d)
        exact=4*d+u*(-4*tr+2*(M+m)*(M*t+m*(1-t)))+4*u*u*D
        out.append(abs(actual-exact)/max(M*M,abs(exact)))
    return max(out)


def vertex(M,m,q2,photon=0.,scale=None):
    scale=M if scale is None else scale
    d=(M*M+m*m-q2)/2
    # The scalar vertex is projected between the on-shell massive spinors.
    def stage(t):
        A=(M*t)**2+(m*(1-t))**2+2*t*(1-t)*d
        L=t*(M*M+d)+(1-t)*(m*m+d)
        N1=-4*L+2*(M+m)*(M*t+m*(1-t))
        N2=4*A
        if photon==0:
            # IR subtraction at the common scale M. O(epsilon) gamma
            # terms multiplying the UV pole are q²-independent.
            J0=4*d/A*np.log(np.sqrt(A)/scale)
            J1=N1/A+N2/(2*A)
            J2=.5*np.log(A/scale**2)-.5
        else:
            l=photon**2
            w=np.sqrt(4*A*l-l*l)
            if w<=0:
                raise ValueError('Photon regulator larger than the test range.')
            U0=np.log(A/l)/(2*A)+(l/(A*w))*(np.arctan2(2*A-l,w)-np.arctan2(-l,w))
            J0=4*d*(U0-np.log(scale/photon)/A)
            def D(u):
                return u*u*A+(1-u)*l
            cuts=np.array([.2,1,4,16,64])*photon/np.sqrt(A)
            cuts=cuts[cuts<1]
            J1=quad(lambda u:(N1+N2*u)*u*u/D(u),0,1,points=cuts,epsabs=3e-13)[0]
            J2=quad(lambda u:u*np.log(D(u)/scale**2),0,1,points=cuts,epsabs=3e-13)[0]
        return -J0-J1-8*J2,4*d/A
    return [quad(lambda t:stage(t)[i],0,1,epsabs=1e-11,epsrel=5e-12)[0] for i in (0,1)]


def radiation(M,m,q2):
    c=(M*M+m*m-q2)/(2*M*m)
    beta=np.sqrt((c-1)*(c+1))/c
    B=2*(np.arctanh(beta)/beta-1)
    def angular(v):
        a=beta*v
        L=np.arctanh(a)/a if a else 1.
        return 2*L-1-(1-beta*beta)/(1-a*a)
    def conversion_integrand(v):
        return (v*angular(v)-B)*v/(1-v*v)
    C=quad(conversion_integrand,0.,1.,epsabs=2e-12,epsrel=3e-12)[0]
    return B,C,conversion_integrand


def vacuum_recoil(M,m,q2,rmax):
    x,w,p,E,B=setup(M,m,np.sqrt(q2),52)
    norm=p/M*4*M*(E+m)
    def G(r):
        return (w@angular_scaled(r,x,M,m,np.sqrt(q2),np.sqrt))/norm
    u0=G(0.)
    return quad(lambda r:(G(r)-u0)/r,0.,rmax,epsabs=1e-12)[0]


def direct_cap_difference(M,m,q2,a,b):
    # Full exact massless real branch in a finite nonzero interval:
    # no logarithm subtraction, vertex or regulator needed.
    def band(s):
        x,w,p,E,B=setup(M,m,np.sqrt(s),64)
        norm=p/M*4*M*(E+m)
        return quad(lambda r:(w@angular_scaled(r,x,M,m,np.sqrt(s),np.sqrt))/(norm*r),
                    a,b,epsabs=2e-11,epsrel=4e-12)[0]
    return band(q2)-band(0.)


def relative_matched_shape(M,m,x,sub,rmax):
    B,C,U=radiation(M,m,x)
    B0,C0,U0=radiation(M,m,0.)
    dv=vacuum_recoil(M,m,x,rmax)-vacuum_recoil(M,m,0.,rmax)
    asym=sub/2+(B-B0)*np.log(rmax)+C-C0+dv
    values=[]
    for lm in (1e-5,1e-7):
        q=M*lm
        v1,w=vertex(M,m,x,q)
        v0,w0=vertex(M,m,0.,q)
        if abs((w-w0)/2-(B-B0))>2e-12:
            raise ValueError('Vacuum real and virtual soft logarithms did not match.')
        v=(v1-v0)/2-(w-w0)/2*np.log(1/lm)
        vmax=np.sqrt((rmax-lm)*(rmax+lm))/rmax
        real=(B-B0)*np.log(rmax/lm)+quad(lambda r:U(r)-U0(r),0.,vmax,epsabs=2e-12)[0]+dv
        values.append(dict(photon_mass_over_M=lm,relative_shape=float(v+real),
                           regulator_defect=float(abs((v+real)-asym))))
    if values[-1]['regulator_defect']>1.6e-6:
        raise ValueError('Finite matched vacuum shape did not approach the regulator-free result.')
    return dict(parent_rest_photon_cap_over_M=rmax,
                shapes_exchange_only_and_real=[float(sub/2),float((B-B0)*np.log(rmax)),float(C-C0),float(dv)],
                inclusive_shape_relative_to_q2_zero=float(asym),regulator_ladder=values)


def run():
    out=[]
    for g in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        rows=[]
        f0,w0=vertex(M,m,0.)
        for q2 in (s*s,.25*(M-m)**2,.8*(M-m)**2):
            f,w=vertex(M,m,q2)
            chi=np.arccosh((M*M+m*m-q2)/(2*M*m))
            ratio=chi/np.tanh(chi)
            algebra=contraction_check(M,m,q2)
            if abs(w/4-ratio)>2e-12 or algebra>2e-13:
                raise ValueError('Either the soft coefficient or the unequal-mass numerator failed.')
            values=[]
            for a in (1e-3,1e-5,1e-7):
                f1,_=vertex(M,m,q2,a*M)
                f2,_=vertex(M,m,0.,a*M)
                values.append(dict(photon_mass_over_M=a,vertex_shape_difference=f1-f2,
                                   cutoff_defect=abs((f1-f2)-(f-f0))))
            if values[-1]['cutoff_defect']>3e-5:
                raise ValueError('Regulator-free relative vertex shape did not reproduce.')
            rows.append(dict(scalar_q2=q2,soft_coefficient=w,soft_coefficient_difference=w-w0,
                             independent_spin_trace_relative_defect=float(algebra),
                             subtracted_vertex_shape_difference=f-f0,
                             photon_mass_ladder=values))
        match=[relative_matched_shape(M,m,s*s,rows[0]['subtracted_vertex_shape_difference'],a)
               for a in (.001,.01,.05)]
        vac_checks=[]
        for a,b in zip(match[:-1],match[1:],strict=True):
            low,high=(v['parent_rest_photon_cap_over_M'] for v in (a,b))
            expected=b['inclusive_shape_relative_to_q2_zero']-a['inclusive_shape_relative_to_q2_zero']
            band=direct_cap_difference(M,m,s*s,low,high)
            discrepancy=abs(expected-band)
            if discrepancy>2e-11:
                raise ValueError('Changing the physical real cut did not reproduce the full real-band test.')
            vac_checks.append(dict(caps=[low,high],absolute_defect=float(discrepancy)))
        backwards=vertex(m,M,s*s,scale=M)
        sym=max(abs(np.asarray(backwards)-vertex(M,m,s*s)))
        if sym>2e-12:
            raise ValueError('The normalized vertex did not exchange its mass labels.')
        out.append(dict(grid=g,rows=rows,vacuum_cap_shape=match,
                        interchange_mass_defect=float(sym),real_band_checks=vac_checks))
    return dict(scope='Feynman-gauge vacuum scalar current shape relative to q²=0, common IR subtraction at M; includes relative vacuum real shape in a declared window, not an absolute or thermal correction.',
                stripped_amplitude_prefactor='g^2/(16 pi^2)',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
