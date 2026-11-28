"""Reduced conservative signed extents for nonparallel shared-medium brackets."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from two_mass_discrete_entropy_gradient_reference import averaged_alphas
from two_mass_entropy_extent_reference import entropy
from two_mass_kinetic_offset_deposition_reference import DT,event
from two_mass_shared_medium_weight_certificate_reference import certificate


def fixture():
    e,m,b,P,g=map(DT,(.046,.002,.007,.36,.012))
    ep=np.hypot(e+m+b,P)
    es,fd=np.asarray([event(e,m,b,P,u)[1]+np.r_[m,b] for u in (-.38,-.31)]).T
    bare=np.r_[ep,min(fd)*DT(.94),fd.mean(),max(fd)*DT(1.06),es]
    B=np.zeros((2,6),dtype=DT)
    B[0,0]=DT(np.pi)**2/(2*ep)
    B[1,1:4]=DT(np.pi)**2/(2*bare[1:4])
    gam=g*g/DT(np.pi)**4
    H=gam*B.T@B
    cap=np.asarray([.9,.8,.85,.9,1.,1.],dtype=DT)
    n=cap*np.asarray([.16,.29,.30,.31,.005,1.1],dtype=DT)
    gn=bare+gam*B.T@(DT(.27)+np.array([.62,.32],dtype=DT)+B@n)
    A=np.array([[-1.,0.,1.,0.,1.,0.],[-1.,1.,0.,0.,0.,1.]],dtype=DT)
    D=np.array([[0.,0.,-1.,1.,0.,0.],[0.,-1.,1.,0.,0.,0.]],dtype=DT)
    return n,cap,gn,H,A,D


def evaluate(xi,n,cap,gn,H,A,D):
    xi=np.asarray(xi,dtype=DT)
    cert=certificate(gn,H,A,D,np.zeros(2,dtype=DT),np.ones(2,dtype=DT),xi,cap,n,np.arange(6)<4)
    if min(cert['map_margin'],cert['occupancy_margin'])<=0 or max(cert['lipschitz'])>=1:
        raise ValueError('The shared fixed-extent control is outside its certified box.')
    w=np.ones(2,dtype=DT)/2
    for i in range(25):
        V=A+w[:,None]*D
        gm=gn+H@(V.T@xi)/2
        ww=-(A@gm)/(D@gm)
        if max(abs(ww-w))<1e-17:
            w=ww
            break
        w=ww
    else:
        raise ValueError('Picard failed.')
    V=A+w[:,None]*D
    dn=V.T@xi
    gm=gn+H@dn/2
    if min(w)<=0 or max(w)>=1 or max(abs(V@gm))>2e-18:
        raise ValueError('Picard failed the individual bracket energy residual.')
    GW=np.diag(D@gm)+(V@H@D.T)*(xi[None,:]/2)
    GX=V@H@V.T/2
    # Solve 2x2 in long double for finite difference near cancellation.
    det=GW[0,0]*GW[1,1]-GW[0,1]*GW[1,0]
    inverse=np.array([[GW[1,1],-GW[0,1]],[-GW[1,0],GW[0,0]]],dtype=DT)/det
    dW=-inverse@GX
    dn_x=V.T+(D.T*xi)@dW
    return w,dW,dn,dn_x,float(max(cert['lipschitz'])),float(max(abs((gn+H@dn)@dn_x)))


n,cap,g0,H,A,D=fixture()
bose=np.arange(6)>=4


def residual(x,h,initial=n):
    gn=g0+H@(initial-n)
    w,dW,dn,dn_x,q,tan=evaluate(x,initial,cap,gn,H,A,D)
    V=A+w[:,None]*D
    bar,He=averaged_alphas(initial,initial+dn,cap,bose)
    db=-He[:,None]*dn_x
    raw=[]; grad=[]; prod=[]
    for j in range(2):
        z=bar[1:4]@V[j,1:4]
        dz=V[j,1:4]@db[1:4]+(bar[1:4]@D[j,1:4])*dW[j]
        y,dy=bar[4+j],db[4+j]
        a=-bar[0]+z+y
        da=-db[0]+dz+dy
        logf=-np.logaddexp(DT(0),bar[0])-np.logaddexp(DT(0),-z)-np.log(-np.expm1(-y))
        dlogf=-db[0]/(1+np.exp(-bar[0]))+dz/(1+np.exp(z))-(np.exp(-y)/(-np.expm1(-y)))*dy
        s=np.sign(a)*np.exp(max(logf,logf-a))*(-np.expm1(-abs(a)))
        raw.append(s);grad.append(s*dlogf+np.exp(logf-a)*da);prod.append(a*s)
    return np.asarray(x,dtype=DT)-h*np.asarray(raw), np.eye(2,dtype=DT)-h*np.asarray(grad),dn,prod,w


def equations(X,h):
    X=np.asarray(X,dtype=DT)
    xi,ww=X[:2],X[2:]
    V=A+ww[:,None]*D
    delta=V.T@xi
    Y=n+delta
    if min(Y)<=0 or max((Y/cap)[~bose])>=1:
        raise ValueError('Full root outside admissible occupancies.')
    a,_=averaged_alphas(n,Y,cap,bose)
    eta=V[:,1:4]@a[1:4]
    z=a[4:]
    # Direct event forward minus reverse, coded separately from the derivative.
    blf=1/(1+np.exp(a[0]))
    bld=1/(1+np.exp(-eta))
    eb=1/(-np.expm1(-z))
    s=blf*bld*eb-(1-blf)*(1-bld)*(eb-1)
    new_gradient=g0+H@delta/2
    return np.asarray(np.r_[xi-h*s,(V@new_gradient)/(D@new_gradient)],float)


def instant(t,x):
    x=np.asarray(x,dtype=DT)
    w=evaluate(np.zeros(2),x,cap,g0+H@(x-n),H,A,D)[0]
    s=-residual(np.zeros(2),DT(1),x)[0]
    return np.asarray((A+w[:,None]*D).T@s,float)
def reduced_step(h,start):
    def f(z): return residual(z,h,start)
    x=root(lambda z:np.asarray(f(z)[0],float),np.zeros(2),
           jac=lambda z:np.asarray(f(z)[1],float),tol=2e-12)
    R,J,delta,prod,w=f(x.x)
    if max(abs(R))>2e-15 or min(prod)<0:
        raise ValueError('Reduced simultaneous step or entropy sign failed.')
    number=max(abs(delta[0]+delta[1:4].sum()),abs(delta[0]+delta[4:].sum()))
    ds=abs(entropy(start+delta,cap,bose)-entropy(start,cap,bose)-h*sum(prod))
    if max(ds,number)>2e-16:
        raise ValueError('Reduced step failed entropy or number conservation.')
    return x.x,J,delta,w,float(ds)


def counted_tangent_test():
    rows=[]
    for v in ((0.,0.),(.003,-.001),(.02,-.01),(-.0025,.007)):
        x=np.asarray(v,dtype=DT)
        out=evaluate(x,n,cap,g0,H,A,D)
        h=DT(1e-5)
        diffs=[[],[]]
        for e in np.eye(2,dtype=DT):
            points=[evaluate(x+j*h*e,n,cap,g0,H,A,D) for j in (-2,-1,1,2)]
            def dv(k):
                return (points[0][k]-8*points[1][k]+8*points[2][k]-points[3][k])/(12*h)
            for row,k in zip(diffs,(0,2),strict=True):
                row.append(dv(k))
        dw,dn=[float(np.max(abs(np.asarray(v).T-ref))) for v,ref in zip(diffs,(out[1],out[3]),strict=True)]
        if dw>2e-12 or dn>5e-14 or out[5]>2e-17:
            raise ValueError('Fixed-extent Jacobian does not reproduce the compensated tangent.')
        rows.append(dict(extents=x.astype(float).tolist(),weights=out[0].astype(float).tolist(),
                         weight_derivative_error=dw,count_derivative_error=dn,
                         contraction_bound=out[4],energy_tangent_defect=out[5]))
    return rows


def steps():
    rows=[]
    for h in map(DT,(.03,.1,.25)):
        x,J,delta,w,ds=reduced_step(h,n)
        # An independently formulated four-unknown problem, using explicit
        # forward and backward factors and no nested weight iteration.
        other=root(equations,np.array([0,0,.27,.68]),args=(h,),tol=1e-12)
        ref=other.x
        dx=max(abs(ref[:2]-x))
        dw=max(abs(ref[2:]-w))
        if max(abs(equations(ref,h)))>1e-14 or max(dx,dw)>1e-14:
            raise ValueError('The reduced and fully coupled equations disagree.')
        xh=np.asarray(x,dtype=DT)
        e=DT(1e-6)
        cols=[]
        for u in np.eye(2,dtype=DT):
            v=[residual(xh+k*e*u,h)[0] for k in (-2,-1,1,2)]
            cols.append((v[0]-8*v[1]+8*v[2]-v[3])/(12*e))
        ej=np.max(abs(J-np.asarray(cols).T))
        de=abs(delta@(g0+H@delta/2))
        if ej>5e-11 or de>2e-17:
            raise ValueError('Reduced extent derivative or energy ledger failed.')
        rows.append(dict(time=float(h),extents=x.tolist(),weights=w.astype(float).tolist(),
                         full_solve_extent_difference=float(dx),full_solve_weight_difference=float(dw),
                         reduced_jacobian_defect=float(ej),energy_identity_defect=float(de),entropy_identity_defect=ds))
    return rows


def trajectory():
    time=.5
    ref=solve_ivp(instant,(0,time),np.asarray(n,float),method='DOP853',
                  rtol=6e-14,atol=3e-16,max_step=.006)
    if not ref.success:
        raise ValueError('High-order trajectory did not reach the benchmark time.')
    rows=[]
    for N in (4,8,16,32):
        u=n.copy()
        ds=0.
        for _ in range(N):
            _,_,du,_,ed=reduced_step(DT(time/N),u)
            u+=du
            ds=max(ds,ed)
        d=u-n
        E=abs(d@(g0+H@d/2))
        rows.append(dict(steps=N,full_energy_defect=float(E),maximum_step_entropy_defect=ds,
                         maximum_count_error=float(max(abs(u-ref.y[:,-1])))))
    order=np.log2(rows[-2]['maximum_count_error']/rows[-1]['maximum_count_error'])
    if abs(order-2)>.06 or max(r['full_energy_defect'] for r in rows)>2e-17:
        raise ValueError('Convergence or global energy ledger failed.')
    return dict(duration=time,order=float(order),mesh=rows)


if __name__=='__main__':
    print(json.dumps(dict(scope='Overlapping but nonparallel brackets in a prescribed reciprocal quadratic medium.',
                          fixed_extent_tangents=counted_tangent_test(),steps=steps(),time_mesh=trajectory()),
                     indent=2,allow_nan=False))
