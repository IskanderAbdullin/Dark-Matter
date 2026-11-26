"""Initial local l=2 pressure response of lowered Bose distributions.
Standalone quadrature, dimensionless normalization of paper 23.
"""
import argparse,json
import numpy as np
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

def occ(E,psi,model):
    s=np.maximum(psi-E,0); b=1/np.expm1(4-psi+E); c=1/np.expm1(4); d=c*(1+c)
    f=np.where(E<psi,b-c,0); fp=np.where(E<psi,b*(1+b),0)
    if model=='quadratic':
        f=np.where(s<1e-5,d*(1+2*c)*s*s/2,f-d*s)
        fp=np.where(s<1e-5,d*(1+2*c)*s,fp-d*(E<psi))
    return f,fp

def response(psi,n,model):
    # Classical dilute MB control with Amax=psi and no individual cutoff.
    mb=model=='MB'; gap=1 if mb else 4-psi
    x,w=leggauss(n); x=(x+1)/2; w=w/2
    la=np.log1p(psi/gap); A=gap*np.expm1(la*x[:,None,None,None]); wa=w[:,None,None,None]*(A+gap)*la
    t=np.pi/4*x[None,:,None,None]; wt=np.pi/4*w[None,:,None,None]
    B=A*np.sin(2*t); uc=np.ones_like(B) if mb else np.clip((psi-A)/B,0,1)
    du=gap+A-B*uc; lu=np.log1p(B*uc/du); u=uc-du*np.expm1(lu*x[None,None,:,None])/B
    wu=w[None,None,:,None]*(du+B*(uc-u))*lu/B
    fun=(lambda E:(np.exp(-E),np.exp(-E))) if mb else lambda E:occ(E,psi,model)
    f1,g1=fun(A+B*u); f2,g2=fun(A-B*u)
    result=0.; signabs=0.
    for lo,hi in [(u,uc),(uc,np.ones_like(uc))]:
        width=np.maximum(0,hi-lo); dv=gap+A-B*hi; lv=np.log1p(B*width/dv)
        v=hi-dv*np.expm1(lv*x[None,None,None,:])/B; wv=w[None,None,None,:]*(dv+B*(hi-v))*lv/B
        f3,g3=fun(A+B*v); f4,g4=fun(A-B*v)
        if mb: d1,d2,d3,d4=f2,f1,-f4,-f3
        else: d1=f2*(1+f3+f4)-f3*f4; d2=f1*(1+f3+f4)-f3*f4; d3=f1*f2-f4*(1+f1+f2); d4=f1*f2-f3*(1+f1+f2)
        for theta in (t,np.pi/2-t):
            P=2*np.sqrt(A)*np.cos(theta); R=2*np.sqrt(A)*np.sin(theta)
            geom=4*A**2.5*np.cos(theta)**2*np.sin(theta)**3*wa*wt*wu*wv
            cov2=(u*v)**2+(1-u*u)*(1-v*v)/2
            def D(sign,ui,vi):
                return ((P*vi/2)**2-(P*ui/2)**2 + sign*P*R/2*(vi*u*v-ui)+R*R/4*(cov2-1))
            D1,D2=D(1,u,v),D(-1,u,v)
            D3,D4=-D(1,v,u),-D(-1,v,u)
            val=geom*(-3*R*R/10)*(g1*d1*D1+g2*d2*D2+g3*d3*D3+g4*d4*D4)
            result+=float(np.sum(val)); signabs+=float(np.sum(abs(val)))
    if mb: pressure=quad(lambda q:q**4*np.exp(-q*q),0,np.inf,epsabs=1e-12)[0]
    else: pressure=quad(lambda q:q**4*occ(q*q,psi,model)[0],0,np.sqrt(psi),epsabs=1e-12)[0]
    return dict(psi=psi,model=model,n=n,pressure=pressure,pressure_derivative=result,rate=result/pressure,absolute_integrand=signabs)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=24);p.add_argument('--model',default='linear');p.add_argument('--psi',type=float,default=3.9);a=p.parse_args(); print(json.dumps(response(a.psi,a.n,a.model)))
