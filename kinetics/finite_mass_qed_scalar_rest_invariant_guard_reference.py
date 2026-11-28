"""Check finite invariant charged regulator in two Lorentz frames."""
import json
from math import fsum

import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_complex_guard_reference import chart

D=np.longdouble
CD=np.clongdouble


def polar(a,b,mid,width,pm):
    h=CD(pm/width)
    z=np.sqrt(1+h*h)
    A,B=h+z,h-z
    if abs(A)>abs(B):B=-1/A
    else:A=-1/B
    l,r=np.exp(np.arcsinh(np.array([a-mid,b-mid],dtype=D)/width))
    H=2*((np.log(r-A)-np.log(l-A))-(np.log(r-B)-np.log(l-B)))/(width*(A-B))
    return H


def audit(M,m,s,cut,t):
    Q0=D(Q)
    M,m,s=map(D,(M,m,s))
    Es,lo,hi,*_=intervals(float(M),float(s),float(m),Q)
    E=D(lo+t*(hi-lo))
    es=np.hypot(s,Q0)
    gam=es/s
    v=Q0/s
    _,_,dat,C,Z,p,w,mass,A,B=chart(M,m,s,E,cut)
    p,w,mass=map(D,(p,w,mass))
    sig=1 if cut=='parent' else -1
    X=(M*M-m*m+sig*s*s)/2
    pstar=np.sqrt((X-mass*s)*(X+mass*s))/s
    e0=X/s
    K=4-2*(M+m)*mass/(Z+M*m)
    wl,wr=e0-D('.0047'),e0+D('.0063')
    pl,pr=(np.sqrt((u-mass)*(u+mass)) for u in (wl,wr))
    def boosted(e):
        r=np.sqrt((e-mass)*(e+mass))
        return gam*e-v*r-w,gam*e+v*r-w
    ends=[*boosted(wl),*boosted(wr),*boosted(e0)]
    if not min(ends)<0<max(ends):raise ArithmeticError('Lab range lost coincidence.')
    peak=[]
    for lim in (wl,wr):
        def hit(q):
            mid,_,_,_,_,_,r=dat(D(q))
            return float(gam*(w+D(q))-v*r*mid-lim)
        peak.append(brentq(hit,float(min(ends)),float(max(ends)),xtol=5e-17))
    ends+=peak
    meta=[]
    for eta in map(D,('1e-6','1e-8','1e-10')):
        eta*=sig
        def s_integrand(delta):
            e=e0+delta
            r=np.sqrt((e-mass)*(e+mass))
            R=pstar*r
            base=e0*e-mass*mass+R
            gap=mass*mass*delta*delta/base
            J=-np.log(base/gap)/(2*R)
            g=-2*sig*s*delta
            return ((2*C+4*Z*J)/(g+1j*eta)+K*J/2*g/(g+1j*eta))*r/2
        def lab_integrand(delta):
            mid,W,um,u0,g,R,r=dat(delta)
            e=w+delta
            lower=max(-1.,(gam*e-wr)/(v*r))
            upper=min(1.,(gam*e-wl)/(v*r))
            if upper<=lower:return 0j
            zr=CD(u0)-1j*eta/g
            H=polar(lower,upper,mid,W,CD(um)-1j*eta/g)
            J=-(np.arcsinh((upper-mid)/W)-np.arcsinh((lower-mid)/W))/(2*R)
            p=np.log(CD(upper)-zr)-np.log(CD(lower)-zr)
            return ((C*p-2*Z*H/R)/g+K*J/2+K*1j*eta*H/(4*R*g))*r/2

        def quad(f,low,high,knots,n):
            X,W=roots_legendre(n)
            a=np.pi*(X+1)/2
            z=np.sin(a/2)**2
            y=np.sin(a)*W*np.pi/4
            pieces=sorted({low,high,*[p for p in knots if low<p<high]})
            raw=[]
            for b,c in zip(pieces[:-1],pieces[1:],strict=True):
                axis=D(b)+(c-b)*z.astype(D)
                raw.append((c-b)*np.dot(y,np.array([f(q) for q in axis])))
            return np.sum(raw)
        knot1=[0.,*[sign*p for p in np.geomspace(1e-16,.006,23) for sign in (-1,1)]]
        knot2=[*ends,0.,*[sign*p for p in np.geomspace(1e-18,.02,87) for sign in (-1,1)],
                *[u+j for u in ends for j in (-.003,-.001,-.0001,-1e-5,-1e-6,-1e-7,1e-7,1e-6,1e-5,.0001,.001,.003)],
                *[b+sg*dd for b in peak for dd in np.geomspace(1e-14,1e-3,22) for sg in (-1,1)]]
        a=quad(s_integrand,wl-e0,wr-e0,knot1,68)
        b=quad(lab_integrand,min(ends),max(ends),knot2,96)
        c=quad(lab_integrand,min(ends),max(ends),knot2,144)
        if abs(a-c)>8e-10 or abs(b-c)>8e-10:
            raise ArithmeticError('The two finite-invariant regulator frames did not close.')
        meta.append(dict(guard=float(eta),scalar_frame=[float(a.real),float(a.imag)],
                         lab_frame=[float(c.real),float(c.imag)],inner_mesh=float(abs(b-c)),
                         frame_difference=float(abs(a-c))))
    return dict(cut=cut,parent_lab_E=float(E),results=meta)


def run():
    a=[audit(.05,.003,s,k,u) for s,k,u in ((.000872,'parent',.36),(.000872,'daughter',.36),
                                         (.000867,'parent',.71),(.000867,'daughter',.71))]
    return dict(scope='Two Lorentz coordinates for a manufactured scalar-rest internal-energy strip and invariant angular regulator. Not a physical medium or width.',controls=a)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
