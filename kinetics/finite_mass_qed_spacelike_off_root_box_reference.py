"""Finite spacelike box with analytically certified exclusion of charged poles."""
import json
import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_spacelike_off_root_patch_reference import panel
from finite_mass_qed_timelike_pole_batched_absorption_reference import load


def wide(lo,hi):
    return np.nextafter(lo,-np.inf),np.nextafter(hi,np.inf)


def plus(a,b):return wide(a[0]+b[0],a[1]+b[1])
def minus(a,b):return wide(a[0]-b[1],a[1]-b[0])


def mul(a,b):
    q=[x*y for x in a for y in b]
    return wide(min(q),max(q))


def frac(a,b):
    if b[0]<=0.:
        raise ArithmeticError('Interval denominator not positive')
    return mul(a,wide(1/b[1],1/b[0]))


def sq(a):
    return wide(0. if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2),max(a[0]**2,a[1]**2))


def sqrt(a):
    if a[0]<=0.:
        raise ArithmeticError('Radicand not certified positive')
    return wide(np.sqrt(a[0]),np.sqrt(a[1]))


def pt(a):return (a,a)


def verify(raw,M,m,Es):
    old=np.diff(np.asarray(raw),axis=1)[:,0]
    pending=[(raw,0)];accepted=[];work=0
    while pending:
        B,depth=pending.pop()
        k,u,c,E=B
        z=None
        try:
            w=mul(k,u)
            t=minus(pt(Es),w)
            Vz=minus(mul(pt(Q),c),k)
            Vx=mul(pt(Q),sqrt(minus(pt(1.),sq(c))))
            V2=plus(sq(Vz),sq(Vx))
            V=sqrt(V2)
            s=minus(sq(t),V2)
            bound=plus(minus(sq(pt(M)),sq(pt(m))),s)
            x=frac(minus(mul(t,E),mul(pt(.5),bound)),V)
            r2=minus(minus(sq(E),sq(pt(M))),sq(x))
            r=sqrt(r2)
            X=frac(Vx,V)
            Z=frac(Vz,V)
            shift=mul(mul(r,X),(-1.,1.))
            pz=plus(mul(Z,x),shift)
            D1=plus(minus(sq(w),sq(k)),mul(pt(2.),minus(mul(w,E),mul(k,pz))))
            sum_D=mul(pt(2.),minus(mul(pt(Es),w),mul(pt(Q),mul(k,c))))
            D2=minus(sum_D,D1)
            if s[1]<0 and all(a*b>0 for a,b in (D1,D2)):
                z=min(abs(x) for q in (D1,D2) for x in q)
        except ArithmeticError:
            pass
        if z is not None:
            accepted.append(z)
            continue
        if depth>=22:raise ArithmeticError(('No uniform off-root box proof at exact subdivisions',B))
        axis=np.argmax(np.diff(np.asarray(B),axis=1)[:,0]/old)
        low,high=B[axis];mid=(low+high)/2.
        for sub in ((low,mid),(mid,high)):
            new=list(B);new[axis]=sub
            pending.append((tuple(new),depth+1))
        work+=1
    return min(accepted),len(accepted),work


def run_box():
    k,u,c=(.001265,.001295),(.493,.507),(-.54,-.46)
    out=[]
    for g in ('fine','third'):
        d,ch,Es=load(g)
        M,m=[x.m for x in ch]
        for ewin in ((1.18,1.26),(1.85,2.25)):
            gap,tested,work=verify((k,u,c,ewin),M,m,Es)
            for channel in ('T','L'):
                vals=[]
                uv=np.concatenate([s.v for s in ch])
                u_knots=np.unique(np.r_[u[0],uv[(uv>u[0])&(uv<u[1])],u[1]])
                for n in (4,6):
                    x,w=roots_legendre(n)
                    steps=[[(a+b)/2+(b-a)*x/2,(b-a)*w/2] for a,b in (k,u,c)]
                    u_segments=[[(a+b)/2+(b-a)*x/2,(b-a)*w/2]
                                for a,b in zip(u_knots[:-1],u_knots[1:],strict=True)]
                    steps[1]=[np.concatenate([row[j] for row in u_segments]) for j in (0,1)]
                    norm=1/(32*np.pi**3*Es)
                    accum=np.zeros(2)
                    for K,W in zip(*steps[0],strict=True):
                        for U,X in zip(*steps[1],strict=True):
                            for C,V in zip(*steps[2],strict=True):
                                r=panel(d,ch,K,K*U,C,channel,ewin,8,24)
                                if r['min_proven_denominator']<gap:
                                    raise ArithmeticError('Box proof violated at a quadrature node')
                                measure=K**3*norm*W*X*V/(4*np.pi)
                                flux=r['raw_off_root_estimate_divided_by_n_scalar']*measure
                                accum+=np.array([flux,-flux*K*U])
                    vals.append(accum)
                ratio=max(abs(vals[0]/vals[1]-1))
                if ratio>1e-7 or min(vals[1][0],-vals[1][1])<=0.:
                    raise ArithmeticError(('Bounded continuum rate check',g,ewin,channel,ratio,vals))
                out.append(dict(grid=g,channel=channel,energy_patch=ewin,certified_subboxes=tested,
                                minimum_abs_denominator=gap,relative_mesh_error=ratio,
                                within_box_velocity_split_count=len(u_knots)-2,
                                normalized_offroot_scalar_number_and_photon_heat=vals[1].tolist()))
    return dict(scope='Finite off-root patch only. Leading collisionless source, no matched inclusive rate.',
                k_range=k,u_range=u,cosine_range=c,data=out)


if __name__=='__main__':
    print(json.dumps(run_box(),indent=2,allow_nan=False))
