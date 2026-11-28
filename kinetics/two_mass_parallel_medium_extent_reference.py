"""Eliminate the weights of parallel energy brackets, then solve the extents."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from two_mass_discrete_entropy_gradient_reference import averaged_alphas
from two_mass_entropy_extent_reference import entropy
from two_mass_kinetic_offset_deposition_reference import DT,event


def fixture():
    e,m,b,P,g=map(DT,(.046,.002,.007,.36,.012))
    M=e+m+b
    dd=np.array([event(e,m,b,P,u)[1]+np.r_[m,b] for u in (-.38,-.31)])
    es,fd=dd.T
    ep=np.hypot(M,P)
    bare=np.array([ep,min(fd)*DT(.94),max(fd)*DT(1.06),*es],dtype=DT)
    B=np.zeros((2,5),dtype=DT)
    B[0,0]=DT(np.pi)**2/(2*ep)
    B[1,1:3]=DT(np.pi)**2/(2*bare[1:3])
    gamma=g*g/DT(np.pi)**4
    H=gamma*B.T@B
    n0=np.asarray([.9,.8,.9,1.,1.],dtype=DT)
    initial=n0*np.asarray([.16,.29,.31,.005,1.1],dtype=DT)
    g0=bare+gamma*B.T@(DT(.27)+np.array([.62,.32],dtype=DT)+B@initial)
    A=np.array([[-1.,1.,0.,1.,0.],[-1.,1.,0.,0.,1.]],dtype=DT)
    d=np.array([0.,-1.,1.,0.,0.],dtype=DT)
    bs=np.asarray([False,False,False,True,True])
    return initial,n0,bs,g0,H,A,d


def chart(xi,n0,cap,bs,g0,H,A,d):
    xi=np.asarray(xi,dtype=DT)
    D=A.T@xi
    c=D@g0+(D@H@D)/2
    b=d@(g0+H@D)
    a=d@H@d/2
    det=b*b-4*a*c
    if det<0 or b<=0:
        raise ValueError('Outside the certified manufactured scalar branch.')
    z=-2*c/(b+np.sqrt(det))
    dn=D+z*d
    gm=g0+H@dn/2
    g1=g0+H@dn
    w=-(A@gm)/(d@gm)
    n1=n0+dn
    if min(w)<=0 or max(w)>=1 or min(n1)<=0 or max((n1/cap)[~bs])>=1:
        raise ValueError('Proposed extents are outside the represented chart.')
    if (d@gm<=0 or d@g1<=0 or abs(xi@w-z)>1e-17 or abs(dn@gm)>2e-18):
        raise ValueError('Aggregate root did not realize the individual energy directions.')
    dv=A.T-np.outer(d,(A@g1)/(d@g1))
    V=A+w[:,None]*d
    dw=-(V@H@dv)/(2*(d@gm))
    return n1,V,dv,dw,z


def evaluate(xi,h,n0,cap,bs,g0,H,A,d):
    n1,V,dv,dw,z=chart(xi,n0,cap,bs,g0,H,A,d)
    avg,Hs=averaged_alphas(n0,n1,cap,bs)
    da=-Hs[:,None]*dv
    x=avg[0]
    S=[]
    aff=[]
    deriv=[]
    for j,u in enumerate(V):
        eta=avg[1:3]@u[1:3]
        y=avg[3+j]
        deta=u[1:3]@da[1:3]+(avg[1:3]@d[1:3])*dw[j]
        Aff=-x+eta+y
        dAff=-da[0]+deta+da[3+j]
        F=-np.logaddexp(DT(0),x)-np.logaddexp(DT(0),-eta)-np.log(-np.expm1(-y))
        sf=1/(1+np.exp(-x))
        sz=1/(1+np.exp(eta))
        dF=-sf*da[0]+sz*deta-np.exp(-y)/(-np.expm1(-y))*da[3+j]
        rate=np.sign(Aff)*np.exp(max(F,F-Aff))*(-np.expm1(-abs(Aff)))
        S.append(rate)
        aff.append(Aff)
        deriv.append(rate*dF+np.exp(F-Aff)*dAff)
    S=np.asarray(S,dtype=DT)
    return (np.asarray(xi,dtype=DT)-h*S,np.eye(len(xi),dtype=DT)-h*np.asarray(deriv,dtype=DT),
            n1,V,z,np.asarray(aff,dtype=DT),S,dv)


def convergence(F):
    initial,cap,bs,g_initial,H,A,d=F
    time=.5
    def instant(t,n):
        n=np.asarray(n,dtype=DT)
        grad=g_initial+H@(n-initial)
        _,_,_,V,_,_,S,_=evaluate(np.zeros(2),1.,n,cap,bs,grad,H,A,d)
        return np.asarray(V.T@S,float)
    ref=solve_ivp(instant,(0,time),np.asarray(initial,float),method='DOP853',
                  rtol=6e-14,atol=3e-16,max_step=.006)
    if not ref.success:
        raise ValueError('Independent continuously energy-matched trajectory failed.')
    out=[]
    for N in (4,8,16,32):
        n=initial.copy()
        h=DT(time/N)
        defects=[]
        counts=[]
        for j in range(N):
            grad=g_initial+H@(n-initial)
            def ev(x):return evaluate(x,h,n,cap,bs,grad,H,A,d)
            c=root(lambda x:np.asarray(ev(x)[0],float),np.zeros(2),jac=lambda x:np.asarray(ev(x)[1],float),
                   tol=2e-12)
            f,_,y,_,_,aff,S,_=ev(c.x)
            if max(abs(f))>2e-15 or min(aff*S)<0:
                raise ValueError('Reduced step failed along the time mesh.')
            ds=entropy(y,cap,bs)-entropy(n,cap,bs)
            defects.append(abs(ds-h*(aff@S)))
            counts.append(c.nfev)
            n=y
        dn=n-initial
        de=abs(dn@(g_initial+H@dn/2))
        err=float(max(abs(n-ref.y[:,-1])))
        if de>2e-17 or max(defects)>2e-16:
            raise ValueError('Long reduced trajectory failed energy or entropy balance.')
        out.append(dict(steps=N,maximum_count_error=err,full_energy_defect=float(de),
                        maximum_entropy_identity_defect=float(max(defects)),
                        max_nonlinear_evaluations=max(counts)))
    order=np.log2(out[-2]['maximum_count_error']/out[-1]['maximum_count_error'])
    if abs(order-2)>.06:
        raise ValueError('The medium-weight trajectory did not converge at second order.')
    return dict(duration=time,order_last_two=float(order),steps=out)


def run():
    F=fixture()
    n0,cap,bs,g0,H,A,d=F
    path=Path(__file__).parent/'results/two-mass-shared-medium-event-weights-2026-11-27.json'
    out=[]
    for r in json.loads(path.read_text())['rows']:
        h=DT(r['step'])
        def ev(x):return evaluate(x,h,*F)
        x=root(lambda x:np.asarray(ev(x)[0],float),np.zeros(2),jac=lambda x:np.asarray(ev(x)[1],float),tol=2e-12)
        fr,J,n1,V,z,aff,S,dv=ev(x.x)
        if max(abs(fr))>2e-15 or min(S*aff)<0:
            raise ValueError('Reduced nonlinear solve or entropy sign failed.')
        checks=[]
        for i in range(2):
            direction=np.eye(2,dtype=DT)[i]
            eps=DT(1e-6)
            zz=np.asarray(x.x,dtype=DT)
            vals=[ev(zz+eps*k*direction)[0] for k in (-2,-1,1,2)]
            deriv=(vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*eps)
            checks.append(max(abs(deriv-J[:,i])))
        delta=n1-n0
        tangent=max(abs((g0+H@delta)@dv))
        diff=max(abs(x.x-r['extents']))
        wd=max(abs(V[:,2]-r['weights']))
        Eerr=abs(delta@(g0+H@delta/2))
        number=max(abs(delta[0]+delta[1:3].sum()),abs(delta[0]+delta[3:].sum()))
        Derr=abs(entropy(n1,cap,bs)-entropy(n0,cap,bs)-h*(aff@S))
        if max(checks)>5e-11 or max(diff,wd,tangent,Eerr,Derr,number)>1e-14:
            raise ValueError('Reduced step failed archive and derivative checks.')
        out.append(dict(step=float(h),extents=x.x.tolist(),weight_root=float(z),weights=V[:,2].astype(float).tolist(),
                        max_jacobian_defect=float(max(checks)),tangent_abs_defect=float(tangent),
                        compared_extent_error=float(diff),compared_weight_error=float(wd),
                        energy_defect=float(Eerr),number_defect=float(number),entropy_defect=float(Derr),solver_evaluations=x.nfev))
    return dict(scope='Manufactured shared quadratic medium with parallel affine energy brackets only.',
                rows=out,trajectory=convergence(F))


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
