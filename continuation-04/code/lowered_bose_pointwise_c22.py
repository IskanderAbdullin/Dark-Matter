"""Pointwise isotropic initial C22 gain/loss, independent q1-fixed quadrature."""
import argparse,json,numpy as np
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes,occupation,local
from lowered_bose_exterior_source import R,V,GYR

def c22(H,psi,q1,n=24):
    if q1*q1>=2*psi or psi<=0:return 0.
    x,w=nodes(n);gap=max(H.eta-psi,1e-15);a=q1*q1;f1=float(occupation(H,psi,a));pmax=np.sqrt(2*psi-a);res=0.
    cuts=np.unique(np.r_[0,min(np.sqrt(psi),pmax),pmax])
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        lg=np.log1p((hi*hi-lo*lo)/(gap+lo*lo));b=lo*lo+(gap+lo*lo)*np.expm1(lg*x[:,None,None]);q2=np.sqrt(b);wb=.5*q2*(b+gap)*lg*w[:,None,None]
        A=(a+b)/2;f2=occupation(H,psi,b)
        c0=np.sqrt(gap*(a+b)/np.maximum(a*b,1e-300));asinh=np.arcsinh(1/c0);c=c0*np.sinh(asinh*x[None,:,None]);wc=c0*np.cosh(asinh*x[None,:,None])*asinh*w[None,:,None]
        B=np.sqrt(np.maximum(0,A*A-a*b*c*c));amin=a*b*c*c/(A+B);rel=.5*(np.sqrt(np.maximum(0,a+b-2*q1*q2*c))+np.sqrt(a+b+2*q1*q2*c))
        uc=np.clip((psi-A)/B,0,1)
        for l,h in [(0,uc),(uc,1)]:
            width=h-l;dv=gap+amin+B*(1-h);lv=np.log1p(B*width/dv);du=dv*np.expm1(lv*x[None,None,:]);q4=amin+B*(1-h)+du;q3=2*A-q4;wv=(dv+du)*lv/B*w[None,None,:]
            f3=occupation(H,psi,q3);f4=occupation(H,psi,q4)
            net=f3*f4*(1+f1+f2)-f1*f2*(1+f3+f4)
            res+=np.sum(wb*wc*wv*rel*net)
    return float(res)

def moments(H,psi,n=24,nq=60):
    x,w=nodes(nq);gap=max(H.eta-psi,1e-15);m=np.zeros(2);abs0=0.;out=0.
    for lo,hi in [(0,psi),(psi,2*psi)]:
        lg=np.log1p((hi-lo)/(gap+lo));q2=lo+(gap+lo)*np.expm1(lg*x);q=np.sqrt(q2);ww=.5*q*(gap+q2)*lg*w
        c=np.array([c22(H,psi,qi,n) for qi in q]);m+=np.array([ww@c,ww@(q2*c)]);abs0+=ww@abs(c)
        if lo>=psi:out=ww@c
    ed=np.array([-psi-1e-12,0,psi+1e-12]);cc,mm,aw=local(H,psi,ed,n)
    return dict(psi=psi,n=n,nq=nq,number=float(m[0]),kinetic=float(m[1]),abs=float(abs0),escape_direct=float(out),escape_signed_half=float(cc[0]/2),escape_ratio=float(out/(cc[0]/2)))

def run(H,n):
    rho=V**2/(4*np.pi*6.6743e-11*R**2);pref=rho*.125*np.sqrt(2)*V*GYR/H.Ic
    results=[]
    for gap in [.5,.1,.01,.001,1e-4,1e-6]:
        psi=H.eta-gap
        for fac in [0,.5,1,3,10]:
            q=np.sqrt(gap)*fac
            if q*q>=psi:continue
            c=c22(H,psi,q,n);f=float(occupation(H,psi,q*q));results.append(dict(gap=gap,q=q,ratio_Gyr=pref*c/f))
    return dict(model=H.model,n=n,pref=pref,rows=results,moments=[moments(H,p,n) for p in [1,3,3.9]])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=24);a=p.parse_args();print(json.dumps(run(Halo(model=a.model),a.n)))
