"""Positive two-moment hard photon closure for a later reciprocal pair step."""
import json
import numpy as np
from scipy.optimize import root
from scipy.special import roots_legendre


class Cell:
    def __init__(self,l,h,n=26):
        self.l=l;self.h=h
        x,w=roots_legendre(n)
        self.p=(l+h)/2+(h-l)*x/2
        self.u=x
        self.W=w*(h-l)/2*self.p*self.p/np.pi**2
        self.B=np.stack([np.ones(n),self.u])
        self.I=self.B*self.W
        self.inv=np.pi**2/self.p

    def f(self,a,b):
        return 1/np.expm1(a+b*self.u)

    def moments(self,f):
        return self.I@f

    def fit(self,Y):
        c=self.W.sum()
        assert Y[0]>0. and abs(Y[1])<Y[0]
        def evaluate(z):
            h=np.exp(z[0])
            a=h*np.cosh(z[1]);b=h*np.sinh(z[1])
            f=self.f(a,b)
            J=-(self.I*(f*(1+f)))@self.B.T
            T=np.array([[a,b],[b,a]])
            return (self.I@f-Y)/c,J@T/c,a,b,J,f
        y0=[np.log(np.log1p(c/Y[0])),0.]
        out=root(lambda z:evaluate(z)[0],y0,jac=lambda z:evaluate(z)[1],tol=2e-11)
        r,_,a,b,J,f=evaluate(out.x)
        assert max(abs(r))<2e-11 and abs(b)<a
        fJ=self.W@(f*self.inv)
        # Derivative of the gauge moment through the Bose entropy fit.
        row=-np.sum(self.I*(f*(1+f)*self.inv),axis=1)
        d=row@np.linalg.inv(J)
        entropy=self.W@(np.log1p(f)+f*np.log1p(1/f))
        return a,b,fJ,d,entropy


def one(cell,a,b):
    F=cell.f(a,b)
    m=cell.moments(F)
    a2,b2,J,d,S=cell.fit(m)
    mag=max(abs(a2-a),abs(b2-b))/max(abs(a),1e-15)
    assert mag<3e-8
    v=m*np.array([.42,-.18])
    h=2e-6
    x,y=[cell.fit(m+j*h*v) for j in (-1,1)]
    numeric=(y[2]-x[2])/(2*h)
    sg=(y[4]-x[4])/(2*h)
    ana=np.array([a2,b2])@v
    assert abs(numeric-d@v)<2e-7*abs(d@v)
    assert abs(sg-ana)<2e-6*max(abs(ana),1e-20)
    return dict(edges=[cell.l,cell.h],a=a,b=b,occupation_range=[min(F),max(F)],
                fit_relative=mag,gauge_derivative_relative=numeric/(d@v)-1,
                entropy_derivative_relative=sg/ana-1)


def run():
    rows=[]
    for l,h in ((.5,.8),(.8,1.3),(1.3,2.6),(2.6,4.5)):
        for a in (.1,1.,6.):
            for s in (-.94,-.5,0.,.5,.94):
                rows.append(one(Cell(l,h),a,s*a))
    return dict(scope='Hard fixed cell, numerical schedule of representable Bose functions; '
                      'no guarantee near condensation or a cell boundary mass deposit.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
