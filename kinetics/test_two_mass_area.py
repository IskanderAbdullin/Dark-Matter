"""Independent adaptive quadrature + exchange/limit/conservation checks."""
import json, time, warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning
from two_mass_area import allowed_area


def bounds(E,m,b,K,D):
    p=np.sqrt(max((E-m)*(E+m),0.))
    return (b*b*E*E+D*D)/(K*E+D*p),(K*E+D*p)/(m*m) if m else np.inf


def reference(e0,e1,r0,r1,x,m,b):
    if x<=m+b or e1==e0 or r1==r0:return 0.
    K=(x*x-m*m-b*b)/2;D=np.sqrt((x-m-b)*(x+m+b)*(x-m+b)*(x+m-b))/2
    pts=[p for r in (r0,r1) if r>0 for p in bounds(r,b,m,K,D) if e0<p<e1]
    def fun(E):
        l,h=bounds(E,m,b,K,D)
        return max(0.,min(r1,h)-max(r0,l))
    # Normalize the integral to rectangle area, not a fixed dimensional atol.
    t=[(p-e0)/(e1-e0) for p in pts]
    return quad(lambda u:fun(e0+u*(e1-e0))/(r1-r0),0,1,points=sorted(set(t)),epsabs=2e-11,epsrel=2e-10,limit=100)[0]


def run():
    rng=np.random.default_rng(53807); rows=[]
    for x,m,b in [(0.1,0.,0.),(.1,.003,0.),(.1,0.,.02),(.1,.003,.02),
                  (1.,.2,.4),(1.,.49,.49),(1.,.499999,.5),(.005,4.3e-12,2e-11),(.005,4.3e-12,.00495),(40.,4e-8,4e-4)]:
        e=np.sort(m+np.exp(rng.uniform(-17,3,(450,2))),axis=1);r=np.sort(b+np.exp(rng.uniform(-17,3,(450,2))),axis=1)
        # Also sample narrow log-grid cells close to the kinematic boundary.
        E=np.geomspace(max(m*1e-7,1e-12),30,120); R=np.geomspace(max(b*1e-7,1e-12),30,120)
        E=np.hypot(E,m);R=np.hypot(R,b)
        i=rng.integers(0,119,250);j=rng.integers(0,119,250)
        e=np.vstack([e,np.stack([E[i],E[i+1]],1)]);r=np.vstack([r,np.stack([R[j],R[j+1]],1)])
        K=(x*x-m*m-b*b)/2;D=np.sqrt((x-m-b)*(x+m+b)*(x-m+b)*(x+m-b))/2
        ii=rng.integers(0,119,250); ee=np.stack([E[ii],E[ii+1]],1)
        rr=[];keep=[]
        for z,ev in enumerate(ee):
            if ev[0]==ev[1]:continue
            l,h=bounds(np.mean(ev),m,b,K,D); q=l if rng.random()<.5 or not np.isfinite(h) else h
            if not np.isfinite(q):continue
            j0=np.searchsorted(R,q)-1
            if 0<=j0<119 and R[j0]<R[j0+1]:keep.append(z);rr.append(R[j0:j0+2])
        if keep:e=np.vstack([e,ee[keep]]);r=np.vstack([r,rr])
        rect=np.diff(e)[:,0]*np.diff(r)[:,0]
        valid=rect>0;e=e[valid];r=r[valid];rect=rect[valid]
        a=allowed_area(e[:,0],e[:,1],r[:,0],r[:,1],x,m,b)/rect
        c=allowed_area(r[:,0],r[:,1],e[:,0],e[:,1],x,b,m)/rect
        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter("always", IntegrationWarning)
            ref=np.array([reference(*ee,*rr,x,m,b) for ee,rr in zip(e,r)])
        err=np.abs(a-ref)
        assert np.max(err)<2e-6,(x,m,b,np.max(err),e[np.argmax(err)],r[np.argmax(err)])
        assert np.max(abs(a-c))<2e-6
        rows.append(dict(x=x,m=m,b=b,cells=len(a),nonzero=int(np.count_nonzero(a)),max_ref_area_error=float(err.max()),max_exchange_error=float(np.max(abs(a-c))),reference_warnings=len(warn)))
    assert allowed_area(1.,2.,1.,2.,1.,.5,.5)==0
    # Subdivision additivity, targeted exact conservation of measure.
    maxadd=0.
    for _ in range(600):
        x=1.;m,b=rng.dirichlet([1,1,1])[:2]
        e=np.sort(m+np.exp(rng.uniform(-9,2,2)));r=np.sort(b+np.exp(rng.uniform(-9,2,2)))
        mid=np.mean(e);midr=np.mean(r);a=allowed_area(*e,*r,x,m,b)
        c=sum(allowed_area(i,j,k,l,x,m,b) for i,j in [(e[0],mid),(mid,e[1])] for k,l in [(r[0],midr),(midr,r[1])])
        maxadd=max(maxadd,float(abs(a-c)/((e[1]-e[0])*(r[1]-r[0]))))
    assert maxadd<1e-7,maxadd
    t=time.perf_counter();n=480;edges=np.geomspace(2.5e-10,40,n+1)
    for x in np.geomspace(.005,40,30):
        E=np.hypot(edges,x*1e-9);R=np.hypot(edges,x*.1)
        a=allowed_area(E[:-1,None],E[1:,None],R[None,:-1],R[None,1:],x,x*1e-9,x*.1)
        assert np.all(a>=0)
    return dict(regimes=rows,subdivision_max_error=maxadd,matrices=30,matrix_n=n,seconds=time.perf_counter()-t)

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
