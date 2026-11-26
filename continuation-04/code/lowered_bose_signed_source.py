"""Initial reversible NR U-U signed count and energy source on binding bins.
Each quadrature event deposits all four legs; no time evolution.
"""
import argparse,json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from lowered_antonov_hessian import Halo
from lowered_bose_exterior_source import R,V,GYR

def nodes(n):
    x,w=leggauss(n);return (x+1)/2,w/2

def occupation(H,psi,q2):
    s=np.maximum(0,psi-q2);y=np.maximum(H.eta-psi+q2,1e-300)
    f=np.exp(y)*np.expm1(s)/(np.expm1(y)*np.expm1(H.eta))
    if H.model=='quadratic':f=np.where(s<1e-5,H.d*(1+2*H.c)*s*s/2+H.d*(1+6*H.c+6*H.c**2)*s**3/6,f-H.d*s)
    return f

def local(H,psi,edges,n=16,logmap=True):
    count=np.zeros(len(edges)-1);moment=np.zeros_like(count);absweight=0.
    if psi<=0:return count,moment,absweight
    x,w=nodes(n);gap=max(H.eta-psi,1e-15)
    if logmap:
        lA=np.log1p(psi/gap);A=gap*np.expm1(lA*x[:,None,None,None]);wa=w[:,None,None,None]*(A+gap)*lA
    else:A=psi*x[:,None,None,None];wa=psi*w[:,None,None,None]
    th=np.pi/4*x[None,:,None,None];wt=np.pi/4*w[None,:,None,None]
    B=A*np.sin(2*th);uc=np.minimum(1,(psi-A)/B)
    if logmap:
        du=np.maximum(gap+A-B*uc,1e-15);lu=np.log1p(B*uc/du);u=uc-du*np.expm1(lu*x[None,None,:,None])/B;wu=w[None,None,:,None]*(du+B*(uc-u))*lu/B
    else:u=uc*x[None,None,:,None];wu=uc*w[None,None,:,None]
    q1=A+B*u;q2=A-B*u;f1=occupation(H,psi,q1);f2=occupation(H,psi,q2);Pu=f1*f2;Su=f1+f2
    geom=8*A**2.5*(np.cos(th)**2*np.sin(th)**3+np.sin(th)**2*np.cos(th)**3)*wa*wt*wu
    for low,high in [(u,uc),(uc,np.ones_like(uc))]:
        width=np.maximum(0,high-low)
        if logmap:
            dv=np.maximum(gap+A-B*high,1e-15);lv=np.log1p(B*width/dv);v=high-dv*np.expm1(lv*x[None,None,None,:])/B;wv=w[None,None,None,:]*(dv+B*(high-v))*lv/B
        else:v=low+width*x[None,None,None,:];wv=width*w[None,None,None,:]
        q3=A+B*v;q4=A-B*v;f3=occupation(H,psi,q3);f4=occupation(H,psi,q4);Pv=f3*f4;Sv=f3+f4
        net=geom*wv*(Pu*(1+Sv)-Pv*(1+Su));absweight+=float(np.sum(abs(net)))
        for q,sign in [(q1,-1),(q2,-1),(q3,1),(q4,1)]:
            s=np.broadcast_to(psi-q,net.shape).ravel();val=sign*net.ravel()
            idx=np.searchsorted(edges,s,side='right')-1
            if np.any(idx<0) or np.any(idx>=len(count)):raise RuntimeError('out of grid')
            count+=np.bincount(idx,weights=val,minlength=len(count));moment+=np.bincount(idx,weights=val*s,minlength=len(count))
    return count,moment,absweight

def run(H,n=16,nr=48,logmap=True):
    e0=H.Qt/H.rt;eta=H.eta
    edges=np.unique(np.r_[np.linspace(-eta,-e0,13),np.linspace(-e0,0,13),np.linspace(0,eta-.2,40),eta-np.geomspace(.2,.001,12),eta]);edges[0]-=1e-12;edges[-1]+=1e-12
    x,w=nodes(nr);r=H.rt*x*x;ww=w*2*H.rt*x*r*r;ps=H.psi(r)
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    c=np.zeros(len(edges)-1);m=np.zeros_like(c);absw=0.;local_err=0.
    for p,weight in zip(ps,ww):
        ci,mi,aw=local(H,p,edges,n,logmap);c+=weight*ci;m+=weight*mi;absw+=weight*aw
        if aw:local_err=max(local_err,abs(ci.sum())/aw,abs(mi.sum())/aw)
    c*=pref;m*=pref;absw*=pref
    cutoff=float(c[edges[1:]<=1e-14].sum());escape=float(c[edges[1:]<=-e0+1e-14].sum());deep=float(c[edges[:-1]>=eta-.2-1e-14].sum())
    return dict(model=H.model,n=n,nr=nr,logmap=logmap,edges=edges.tolist(),count=c.tolist(),binding_moment=m.tolist(),absweight=absw,count_ledger=float(c.sum()),binding_ledger=float(m.sum()),local_ledger_max=local_err,cutoff=cutoff,escape=escape,deep_gap_02=deep)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=16);p.add_argument('--nr',type=int,default=48);p.add_argument('--uniform',action='store_true');a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr,not a.uniform)))
