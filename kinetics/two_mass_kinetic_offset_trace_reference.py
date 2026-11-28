"""Independent two-mass geometric energy area and trace moments, in kinetic offsets."""
import argparse
import json

import numpy as np
from scipy.special import roots_legendre

DT=np.longdouble


def offset(t,m):
    return 2*m*np.sinh(t/2)**2


def delta_offset(a,b,m):
    return 2*m*np.sinh((a+b)/2)*np.sinh((b-a)/2)


def log_remainder(z):
    # z-log1p(z): needed when both absolute traces round to the rest mass.
    if abs(z)<1e-3:
        return z*z*(1/2+z*(-1/3+z*(1/4+z*(-1/5+z*(1/6-z/7)))))
    return z-np.log1p(z)


def geometry(parent_excess,m,b,ps,pf,order=32):
    """All outputs are normalized by the *whole* energy-cell rectangle.

    ps and pf are momentum-edge pairs in the same units as the positive
    daughter masses. Specify the parent rest excess directly; a near
    threshold excess lost when rounding the parent mass is unrecoverable.
    """
    if not np.all(np.isfinite([parent_excess,m,b,*ps,*pf])):
        raise ValueError("A finite two-mass rectangle is required.")
    e,m,b=map(DT,(parent_excess,m,b))
    if min(e,m,b)<=0 or min(ps[0],pf[0])<0 or ps[1]<=ps[0] or pf[1]<=pf[0]:
        raise ValueError("Use strictly positive masses/excess and ordered nonnegative momentum edges.")
    Y=2*np.arcsinh(np.sqrt(e*(2*(m+b)+e)/(4*m*b)))
    lo,hi=np.arcsinh(np.asarray(ps,dtype=DT)/m)
    a,z=np.arcsinh(np.asarray(pf,dtype=DT)/b)
    ds,dr=delta_offset(lo,hi,m),delta_offset(a,z,b)
    if min(ds,dr)<=0:
        raise ValueError("Rapidity edge resolution was insufficient.")
    edges=sorted({a,z,*(x for k in (lo,hi) for x in (Y+k,abs(Y-k)) if a<x<z)})
    v,w=(np.array(x,dtype=DT) for x in roots_legendre(order))
    out=np.zeros(5,dtype=DT)
    for left,right in zip(edges[:-1],edges[1:],strict=True):
        # Integrate over dimensionless fermion rapidity. Stable widths
        # remain nonzero when the absolute E edges coincide in float64.
        for node,weight in zip(v,w,strict=True):
            t=left+(right-left)*(1+node)/2
            s0,s1=max(lo,abs(Y-t)),min(hi,Y+t)
            if s1<=s0:
                continue
            width=delta_offset(s0,s1,m)
            kr=offset(t,b)
            ks=offset(s0,m)
            j=b*np.sinh(t)*(right-left)*weight/(2*dr)
            f=width/ds
            g=width/(m+ks)
            scalar_trace=m*m*np.log1p(g)/ds
            scalar_deficit=m*(ks*g+m*log_remainder(g))/ds
            out+=j*np.array([f,scalar_trace,f*b*b/(b+kr),scalar_deficit,f*b*kr/(b+kr)],dtype=DT)
    return np.asarray(out,float)


def mp_reference(e,m,b,ps,pf):
    # Independent high-precision quadrature of the original E±
    # rectangle with K r ± D sqrt(r²-b²), no kinetic-offset primitives.
    import mpmath as mp
    with mp.workdps(85):
        M,B,Q=map(lambda z:mp.mpf(float(z)),(m,b,e))
        ell,hi=[mp.sqrt(M*M+mp.mpf(float(x))**2) for x in ps]
        a,z=[mp.sqrt(B*B+mp.mpf(float(x))**2) for x in pf]
        K=M*B+Q*(2*(M+B)+Q)/2
        D=mp.sqrt((K-M*B)*(K+M*B))
        # Inverse kinematics: r bounds at each fixed scalar energy edge.
        def cut(x):
            p=mp.sqrt((x-M)*(x+M))
            return ((K*x-D*p)/M**2,(K*x+D*p)/M**2)
        cuts=sorted({a,z,*(v for x in (ell,hi) for v in cut(x) if a<v<z)})
        whole=(hi-ell)*(z-a)
        def rec(r,index):
            p=mp.sqrt(max(mp.mpf(0),(r-B)*(r+B)))
            x=max(ell,(K*r-D*p)/B**2)
            y=min(hi,(K*r+D*p)/B**2)
            if y<=x:
                return mp.mpf(0)
            d=y-x
            tr=M*M*mp.log(y/x)
            return (d,tr,B*B*d/r,M*d-tr,B*d-B*B*d/r)[index]
        return np.array([float(mp.quad(lambda x:rec(x,i),cuts)/whole) for i in range(5)])


def width(p,m):
    a,h=np.arcsinh(np.array(p,dtype=DT)/DT(m))
    return float(delta_offset(a,h,DT(m)))


def subdivision(e,m,b,ps,pf,whole):
    x=(ps[0]+ps[1])/2
    y=(pf[0]+pf[1])/2
    L=width(ps,m)*width(pf,b)
    value=np.zeros(5)
    for a in ((ps[0],x),(x,ps[1])):
        for c in ((pf[0],y),(y,pf[1])):
            value+=geometry(e,m,b,a,c,64)*(width(a,m)*width(c,b)/L)
    return float(np.max(abs(whole-value)/np.maximum(abs(whole),1e-100)))


def tests():
    configs=[
      ("ordinary_partial",.093,.061,.122,[.024,.081],[.02,.15]),
      ("near_threshold_double_collapse",2e-23,.061,.122,[1e-12,5e-12],[.4e-12,7e-12]),
      ("asymmetric_one_near_rest",.031,.065,.018,[8e-13,7e-11],None),
      ("full",.031,.061,.122,[.051,.060],[.10,.13]),
      ("empty",.09,.05,.08,[0.,.001],[0.,.001]),
    ]
    cases=[]
    for name,e,m,b,ps,pf in configs:
        if pf is None:
            em,mm,bb=map(DT,(e,m,b))
            Y=2*np.arcsinh(np.sqrt(em*(2*(mm+bb)+em)/(4*mm*bb)))
            p=b*np.sinh(Y)
            pf=[float(p-1e-11),float(p+3e-11)]
        low,high=geometry(e,m,b,ps,pf),geometry(e,m,b,ps,pf,64)
        ref=mp_reference(e,m,b,ps,pf)
        tol=3e-7 if name=="asymmetric_one_near_rest" else 3e-11
        compare=float(np.max(abs(high-ref)/np.maximum(abs(ref),1e-100)))
        if compare>tol or np.max(abs(low-high)/np.maximum(abs(high),1e-100))>tol:
            raise ValueError(f"Offset quadrature or independent reference failed at {name}: {compare}; {low}, {high}, {ref}.")
        closure=abs(high[1]+high[3]-m*high[0])+abs(high[2]+high[4]-b*high[0])
        if not 0<=high[0]<=1+2e-14 or closure>3e-15*max(m,b):
            raise ValueError(f"Trace/deficit closure at {name} failed.")
        other=geometry(e,b,m,pf,ps,64)[[0,2,1,4,3]]
        if np.max(abs(high-other)/np.maximum(abs(high),1e-100))>tol:
            raise ValueError(f"Daughter exchange failed at {name}.")
        collapsed=[np.hypot(w,p[0])==np.hypot(w,p[1]) for w,p in ((m,ps),(b,pf))]
        if name in ("near_threshold_double_collapse","asymmetric_one_near_rest") and (not any(collapsed) or high[0]<=0):
            raise ValueError("The intended zero-width absolute-energy control was not exercised.")
        sub=subdivision(e,m,b,ps,pf,high)
        if sub>tol or (name=='full' and abs(high[0]-1)>1e-14):
            raise ValueError(f"Exact cell split or full cell failed at {name}: {sub}.")
        cases.append(dict(name=name,parent_excess=e,scalar_mass=m,fermion_mass=b,
                          parent_excess_lost_in_absolute_mass=(m+b+e==m+b),
                          kinetic_energy_cell_widths=[width(ps,m),width(pf,b)],
                          momentum_edges=[ps,pf],absolute_energy_width_collapsed=list(map(bool,collapsed)),
                          normalized_geometry_and_traces=high.tolist(),reference=ref.tolist(),
                          relative_error=compare,normalized_trace_defect_closure=float(closure),
                          subdivision_relative_error=sub))
    return cases


def check_upstream():
    from two_mass_area import allowed_area
    rng=np.random.default_rng(7419)
    errors=[]
    count=0
    for i in range(96):
        m,b=np.exp(rng.uniform(-4.,0.,2))
        Y=rng.uniform(.08,3.7)
        M=np.sqrt(m*m+b*b+2*m*b*np.cosh(Y))
        e=M-m-b
        t=rng.uniform(0.,4.,2)
        ds,df=rng.uniform(.01,.6,2)
        ps=m*np.sinh([t[0],t[0]+ds])
        pf=b*np.sinh([t[1],t[1]+df])
        h=np.hypot(m,ps)
        f=np.hypot(b,pf)
        upstream=allowed_area(*h,*f,M,m,b).item()/((h[1]-h[0])*(f[1]-f[0]))
        v=geometry(e,m,b,ps,pf)[0]
        count+=0<v<1.-2e-15
        errors.append(float(abs(v-upstream)))
    if max(errors)>2e-9 or count<20:
        raise ValueError(f"Ordinary-width upstream comparison failed: {max(errors)}, mixed={count}.")
    return dict(controls=len(errors),mixed=int(count),max_absolute_area_fraction_discrepancy=max(errors))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check-upstream',action='store_true')
    a=p.parse_args()
    print(json.dumps(dict(scope="Independent geometric two-mass shell-weight control, not a collision rate or evolution.",
                          fields=["area_fraction","scalar_trace","fermion_trace","scalar_trace_deficit","fermion_trace_deficit"],
                          cases=tests(),ordinary_width_upstream=(check_upstream() if a.check_upstream else None)),
                     indent=2,allow_nan=False))
