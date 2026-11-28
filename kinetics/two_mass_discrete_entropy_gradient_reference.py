"""Second-order unsplit event reference with an exact discrete entropy gradient."""
import json

import numpy as np
from scipy.integrate import solve_ivp

from two_mass_coupled_extent_newton_reference import solve as backward
from two_mass_entropy_extent_reference import DT, batch, entropy, rates


def meanlog(a,b):
    """Mean of log along a positive line segment and its end derivative."""
    v=(b-a)/a
    if abs(v)<DT(1e-3):
        ser=DT(0)
        slope=DT(0)
        for j in range(7,0,-1):
            ser=ser*v+DT((-1)**(j+1))/(j*(j+1))
            slope=slope*v+DT((-1)**(j+1))/(j+1)
        return np.log(a)+v*ser,slope/a
    ratio=b/a
    log=np.log(ratio)
    return np.log(a)+(ratio*log-v)/v,(v-log)/(a*v*v)


def averaged_alphas(before,after,g,bose):
    bar=np.empty_like(before)
    H=np.empty_like(before)
    for i,(x,y,c,is_bose) in enumerate(zip(before,after,g,bose,strict=True)):
        q,dq=meanlog(x,y)
        if is_bose:
            r,dr=meanlog(c+x,c+y)
            H[i]=dq-dr
        else:
            r,dr=meanlog(c-x,c-y)
            H[i]=dq+dr
        bar[i]=r-q
    if min(H)<=0:
        raise ValueError("Invalid negative Hessian for the averaged entropy coordinates.")
    return bar,H


def flux(begin,end,g,bose,V):
    a,H=averaged_alphas(begin,end,g,bose)
    ss=[]
    aa=[]
    ds=[]
    for u in V.T:
        x=a[0]
        y=(u*a)@bose
        z=(u*a)@~bose-u[0]*a[0]
        A=u@a
        unblocked=np.exp(-np.logaddexp(DT(0),-x))
        Ff=np.exp(-np.logaddexp(DT(0),z))
        B=np.exp(-y)/(-np.expm1(-y))
        fwd=-np.logaddexp(DT(0),x)-np.logaddexp(DT(0),-z)-np.log(-np.expm1(-y))
        S=np.sign(A)*np.exp(max(fwd,fwd-A))*(-np.expm1(-abs(A)))
        rev=unblocked*Ff*B
        dx=(B*u*bose-Ff*u*~bose)*H
        dx[0]=unblocked*H[0]
        ds.append(S*dx-rev*u*H)
        aa.append(A)
        ss.append(S)
    return np.array(ss,dtype=DT),np.array(aa,dtype=DT),np.array(ds,dtype=DT)


def solve(begin,g,bose,V,h):
    x=np.zeros(V.shape[1])
    upper=np.where(bose,np.inf,g)
    def valid(y):
        return np.all(y>0) and np.all(y<upper)
    for j in range(30):
        end=begin+V@x
        ss,A,G=flux(begin,end,g,bose,V)
        f=x-h*np.asarray(ss,float)
        err=np.linalg.norm(f,np.inf)
        if err<3e-13:
            change=entropy(end,g,bose)-entropy(begin,g,bose)
            predicted=h*(ss@A)
            if abs(change-predicted)>3e-11 or change<-3e-14 or predicted<0.:
                raise ValueError("Discrete entropy equality or production sign did not close.")
            return end,x,float(abs(change-predicted)),j
        matrix=np.eye(len(x))-h*np.asarray(G@V,float)
        try:
            dx=np.linalg.solve(matrix,-f)
        except np.linalg.LinAlgError:
            return None
        m=1.
        while m>=2.**-30:
            z=x+m*dx
            q=begin+V@z
            if valid(q):
                news=flux(begin,q,g,bose,V)[0]
                if np.linalg.norm(z-h*np.asarray(news,float),np.inf)<err*(1-1e-4*m):
                    x=z
                    break
            m*=.5
        else:
            return None
    return None


def one(title,e,m,b,p,T):
    K,g,bs,V=batch(e,m,b,p,(-.65,.21,.72))
    V=np.asarray(V,dtype=DT).T
    eta=np.where(bs,DT(-.4),DT(-.2))
    eta[0]=DT(-.4)+DT(-.2)-DT(e)/DT(T)
    initial=g/(np.exp(K/DT(T)-eta)+np.where(bs,-1,1))
    initial[0]=DT(.7)*g[0]
    def exact_rhs(t,y):
        return np.array(V@np.array([rates(y,g,bs,v)[1] for v in V.T]),float)
    time=.8
    u=solve_ivp(exact_rhs,(0.,time),np.array(initial,float),method="DOP853",rtol=5e-13,atol=2e-15,max_step=.008)
    if not u.success:
        raise ValueError("Unsplit independent flow failed.")
    ref=u.y[:,-1]
    out=[]
    f=(~bs)&(np.arange(len(K))>0)
    for N in (8,16,32,64):
        n=initial.copy()
        max_entropy=0.
        maxit=0
        for z in range(N):
            one=solve(n,g,bs,V,time/N)
            if one is None:
                raise ValueError("Unexpected nonconvergence in the smooth regime.")
            new,ex,de,it=one
            delta=new-n
            mass=float(abs((K@delta-DT(e)*sum(ex))/max(K)))
            number=max(abs(delta[bs].sum()+delta[0]),abs(delta[f].sum()+delta[0]))
            if max(mass,number)>5e-15 or min(new)<=0 or max((new/g)[~bs])>=1:
                raise ValueError("Energy, number or occupation ledger failed.")
            n=new
            max_entropy=max(max_entropy,de)
            maxit=max(maxit,it)
        out.append(dict(steps=N,maximum_count_error=float(max(abs(n-ref))),
                        maximum_entropy_identity_defect=max_entropy,maximum_newton_iterations=maxit))
    order=np.log2(out[-2]['maximum_count_error']/out[-1]['maximum_count_error'])
    if abs(order-2)>0.05:
        raise ValueError("Second order not reproduced.")
    # A separate guarded severe step; if the averaged-affinity solve
    # refuses the step, use the previously checked backward scheme.
    initial[0]=g[0]*DT(.9999)
    initial[1:]=g[1:]*DT(1e-10)
    big=solve(initial,g,bs,V,500.)
    substeps=[]
    if big is None:
        def fallback(start,h,level):
            try:
                final=backward(start,g,bs,V,h)[0]
            except (ValueError,np.linalg.LinAlgError):
                if level>=10:
                    raise ValueError("No complete feasible integration at fallback limit.")
                return fallback(fallback(start,h/2,level+1),h/2,level+1)
            substeps.append(h)
            return final
        finish=fallback(initial,500.,0)
        used="backward_fallback"
    else:
        finish=big[0]
        used="discrete_entropy_gradient"
        substeps=[500.]
    if abs(sum(substeps)-500.)>1e-13:
        raise ValueError("Fallback did not cover the requested interval.")
    if entropy(finish,g,bs)<entropy(initial,g,bs):
        raise ValueError("Guarded stiff step lowered entropy.")
    return dict(title=title,last_mesh_order=float(order),mesh=out,stiff_solver=used,stiff_substeps=substeps,
                stiff_entropy_increase=float(entropy(finish,g,bs)-entropy(initial,g,bs)))


if __name__=='__main__':
    print(json.dumps(dict(scope="Manufactured fixed-chart entropy gradient and time integration, no physical dynamics.",
        cases=[one('ordinary',.093,.061,.122,.33,.07),one('near-rest',2e-23,.061,.122,4e-12,4e-23)]),
                     indent=2,allow_nan=False))
