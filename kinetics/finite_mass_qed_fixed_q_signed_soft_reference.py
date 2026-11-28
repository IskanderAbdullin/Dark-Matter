"""Fixed-q two-body recoil phase space with leading signed soft subtraction."""
import gzip,json
from pathlib import Path
import bisect
import functools
import mpmath as mp
import numpy as np
from scipy.special import roots_legendre
from early_pole_decay_mode_clock_reference import intervals, Q

base=Path(__file__).parent/'results'

@functools.lru_cache(maxsize=16)
def legendre(n):
    return roots_legendre(n)


def F(M,m,s,ep,q,k,n=64,naz=32):
    M,m,s,ep,q,k=map(np.longdouble,(M,m,s,ep,q,k))
    E=np.hypot(s,q)
    p=np.sqrt((ep-M)*(ep+M))
    ell0=M*M/(ep+p)
    ell1=ep+p
    eb=ep-k-E
    if eb<=m:
        return 0.,0.,[ell0,ell1]
    kb=np.sqrt((eb-m)*(eb+m))
    A0=(ep-k)*E-(M*M+s*s-m*m)/2
    # t=k*(E_p-p z); solve the squared support condition analytically.
    if k!=0:
        b=np.sort((q*q-A0+q*kb*np.array([-1.,1.]))/k)
        l,r=max(ell0,b[0]),min(ell1,b[1])
    else:
        l,r=ell0,ell1
        if abs(A0)>p*q:
            return 0.,0.,[l,r]
    if r<=l:
        return 0.,0.,[l,r]
    x,w=legendre(n)
    shift=m*m/(eb+kb)
    start,end=np.log([l-ell0+shift,r-ell0+shift])
    nodes=(start+end)/2+(end-start)*x/2
    t=np.exp(nodes)-shift
    ell=ell0+t
    weight=w*(end-start)*np.exp(nodes)/(4*p)
    zz=1-t/p
    vv=np.sqrt(np.maximum(0.,(t/p)*(1+zz)))
    q2=np.sqrt((p-k)*(p-k)+2*k*t)
    cc=(A0+k*ell)/(q*q2)
    if np.max(abs(cc))>1.00000001:
        raise ValueError('Collinear analytic interval did not bound the on-shell direction.')
    cc=np.clip(cc,-1.,1.)
    zzQ=(p*zz-k)/q2
    vQ=p*vv/q2
    phi=np.cos((np.arange(naz)+.5)*2*np.pi/naz)
    sp=q*(cc[:,None]*zzQ[:,None]+(np.sqrt((1-cc)*(1+cc))*vQ)[:,None]*phi)
    bn=ell[:,None]-E+sp
    a=-ell[:,None]; b=-bn
    d=(M*M+m*m-s*s)/2+k*(bn-ell[:,None])
    V=(1/a-1/b)/2
    H=M*M/(a*a)+m*m/(b*b)-2*d/(a*b)
    X=b*M*M/a+a*m*m/b-2*d
    amp=-4*(d+M*m)*H-8*k*V*X+16*k*k*V*V*a*b
    if np.min(amp)<-1e-9:
        raise ValueError('Negative real photon norm.')
    norm=2*((M+m)**2-s*s)
    value=weight@(amp.mean(axis=1)*p/(q2*norm))
    return value,(r-l)/(ell1-ell0),[float(l),float(r)]


def chart_b(M,m,s,ep,q,k,n=72,naz=36):
    if k==0.:
        raise ValueError('The inverse eliminated chart requires a nonzero photon.')
    es=np.hypot(s,q)
    P=np.sqrt((ep-M)*(ep+M))
    E=ep-es
    eb=E-k
    kb=np.sqrt((eb-m)*(eb+m))
    L=M*M+s*s-m*m-2*ep*es
    D=k*E-L/2
    v0=D-k*k-k*kb
    bound=np.sort((np.array([-P*q,P*q])-v0)/k)
    l,h=max(0.,bound[0]),min(2*kb,bound[1])
    if h<=l:
        return 0.
    shift=m*m/(eb+kb)
    L,H=np.log([l+shift,h+shift])
    z,w=legendre(n)
    log=(L+H)/2+(H-L)*z/2
    v=np.exp(log)-shift
    wt=w*(H-L)*np.exp(log)/4
    u=v0+k*v
    Q=np.sqrt(P*P+q*q-2*u)
    cs=(k+kb-v)/Q
    if max(abs(cs))>1+1e-9:
        raise ValueError('Inverse elimination selected an invalid emission angle.')
    cs=np.clip(cs,-1.,1.)
    ys=np.clip(u/(P*q),-1.,1.)
    C=(u-q*q)/Q
    S=P*q*np.sqrt((1-ys)*(1+ys))/Q
    ph=np.cos((np.arange(naz)+.5)*2*np.pi/naz)
    sp=(C*cs)[:,None]+(S*np.sqrt((1-cs)*(1+cs)))[:,None]*ph
    Bn=(shift+v)[:,None]
    ell=Bn+es-sp
    a,b=-ell,-Bn
    d=(M*M+m*m-s*s)/2+k*(Bn-ell)
    V=(1/a-1/b)/2
    U=M*M/(a*a)+m*m/(b*b)-2*d/(a*b)
    X=b*M*M/a+a*m*m/b-2*d
    amp=-4*(d+M*m)*U-8*k*V*X+16*k*k*V*V*a*b
    return float(wt@(amp.mean(axis=1)/(Q*2*((M+m)**2-s*s))))


def convergence_coeffs(M,m,s,ep,q,n=44,naz=12):
    # All angles are allowed in an open interval around k=0 for
    # the interior Born parents tested here. Differentiate only there.
    with mp.workdps(61):
        M,m,s,ep,q=map(mp.mpf,(M,m,s,ep,q))
        E=mp.sqrt(s*s+q*q)
        p=mp.sqrt((ep-M)*(ep+M))
        ell0=M*M/(ep+p)
        ell1=ep+p
        xb=ep-E
        shift=m*m/(xb+mp.sqrt((xb-m)*(xb+m)))
        l,h=mp.log(shift),mp.log(ell1-ell0+shift)
        norm=2*((M+m)**2-s*s)
        x,w=legendre(n)
        angles=[]
        for xx,ww in zip(x,w,strict=True):
            a=mp.exp((h+l)/2+(h-l)*mp.mpf(xx)/2)
            t=a-shift
            zz=1-t/p
            vv=mp.sqrt((t/p)*(1+zz))
            angles.append((ell0+t,t,zz,vv,mp.mpf(ww)*(h-l)*a/(4*p)))
        phi=[mp.cos((i+mp.mpf('.5'))*2*mp.pi/naz) for i in range(naz)]
        def func(k):
            A0=(ep-k)*E-(M*M+s*s-m*m)/2
            result=[]
            for ell,t,z,v,weight in angles:
                cq=mp.sqrt((p-k)**2+2*k*t)
                cs=(A0+k*ell)/(q*cq)
                ss=mp.sqrt(1-cs*cs)
                nv=(p*z-k)/cq
                nt=p*v/cq
                tmp=[]
                for cp in phi:
                    bn=ell-E+q*(cs*nv+ss*nt*cp)
                    a,b=-ell,-bn
                    d=(M*M+m*m-s*s)/2+k*(bn-ell)
                    V=(1/a-1/b)/2
                    H=M*M/(a*a)+m*m/(b*b)-2*d/(a*b)
                    X=b*M*M/a+a*m*m/b-2*d
                    tmp.append(-4*(d+M*m)*H-8*k*V*X+16*k*k*V*V*a*b)
                result.append(weight*mp.fsum(tmp)*p/(cq*norm*naz))
            return mp.fsum(result)
        C=[mp.diff(func,mp.mpf(0),j)/mp.factorial(j) for j in range(5)]
        e=mp.mpf('1e-7')
        exact=[func(e)-C[0],func(e)+func(-e)-2*C[0]]
        fitted=[sum(C[j]*e**j for j in range(1,5)),2*(C[2]*e**2+C[4]*e**4)]
        err=max(abs(a/b-1) for a,b in zip(fitted,exact,strict=True))
        if err>mp.mpf('1e-6'):
            raise ValueError('Zero-cut expansion did not reproduce direct high-precision differences.')
        return np.asarray(C,float),float(err)


def photon_boundaries(M,m,s,ep,q):
    E=np.hypot(s,q)
    P=np.sqrt((ep-M)*(ep+M))
    out=[]
    xb=ep-E
    for z in (-1,1):
        for c in (-1,1):
            t=P-c*q
            A=(xb-t)*(xb+t)-m*m
            B=2*(xb-z*t)
            out.append(A/B)
    return np.asarray(out)


def quantum_setup(d,M,m,ep,E):
    eped=np.hypot(d['endpoint_parent_momentum_edges'],M).tolist()
    ebed=np.hypot(d['endpoint_fermion_momentum_edges'],m).tolist()
    fp=d['endpoint_parent_occupations'][bisect.bisect_right(eped,ep)-1]
    def f(k):
        return d['endpoint_fermion_occupations'][bisect.bisect_right(ebed,ep-k-E)-1]
    return fp,f,np.asarray(ebed)


def one(data,M,m,s,ep,q):
    E=np.hypot(s,q)
    y0,_,_=F(M,m,s,ep,q,0.)
    signed=(1e-7,1e-6,3e-6,1e-5)
    chart_error=float(max(abs(F(M,m,s,ep,q,zz)[0]-chart_b(M,m,s,ep,q,zz))
                          for z in signed for zz in (-z,z)))
    if chart_error>3e-12:
        raise ValueError('The two independent energy-delta angle charts disagree.')
    @functools.lru_cache(maxsize=12)
    def zero(n,a): return F(M,m,s,ep,q,0.,n,a)[0]
    chi=np.arccosh((M*M+m*m-s*s)/(2*M*m))
    reference=2*(chi/np.tanh(chi)-1)
    if abs(y0-reference)>4e-12:
        raise ValueError('Fixed-q on-shell normalization failed its eikonal limit.')
    c,cd=convergence_coeffs(M,m,s,ep,q)
    if abs(c[0]-y0)>1e-11:
        raise ValueError('High precision derivative mesh lost the leading eikonal coefficient.')
    fp,fd,b_edges=quantum_setup(data,M,m,ep,E)
    same=fd(0.)
    wp=fp*(1-same)
    wm=(1-fp)*same
    stim=wp-wm
    limits=np.r_[abs(photon_boundaries(M,m,s,ep,q)),abs(ep-E-b_edges)]
    if min(limits)<2e-7:
        raise ValueError('The fixture is too close to a native occupation or kinematic cut.')
    tiny=1e-7

    def integrands(k,n,az):
        nB=1/np.expm1(k/.9)
        if k<=tiny:
            o=c[1]*k+c[3]*k**3
            ev=c[2]*k*k+c[4]*k**4
            return np.array([((wp+wm)*o+(wp-wm)*ev)/k,2*nB*stim*ev/k])
        y0=zero(n,az)
        p,_,_=F(M,m,s,ep,q,k,n,az)
        a,_,_=F(M,m,s,ep,q,-k,n,az)
        bp,bm=fd(k),fd(-k)
        vac=(fp*(1-bp)*(p-y0)-(1-fp)*bm*(a-y0))
        vac+=y0*(fp*(same-bp)+(1-fp)*(same-bm))
        th=stim*(p+a-2*y0)+(same-bp)*p+(same-bm)*a
        return np.array([vac,nB*th])/k

    def integral(lo,hi,n,a,outer):
        t,w=legendre(outer)
        points=sorted({z for z in np.r_[lo,limits,tiny,hi] if lo<=z<=hi})
        val=np.zeros(2,dtype=np.longdouble)
        for u,v in zip(points[:-1],points[1:],strict=True):
            if u>=tiny and v/u>1.2:
                l,r=np.log([u,v])
                x=np.exp((l+r)/2+(r-l)*t/2)
                wt=w*(r-l)*x/2
            else:
                x=(u+v)/2+(v-u)*t/2
                wt=w*(v-u)/2
            val+=wt@np.asarray([integrands(k,n,a) for k in x])
        return np.asarray(val,float)

    mesh=[]
    for L in (2.5e-6,5e-6,1e-5):
        at=integral(0.,L,64,24,80)
        small=integral(0.,L,40,12,40)
        mini=integral(1e-10,L,64,24,80)
        if max(abs(at-small))>2e-8:
            raise ValueError('The signed fixed-q real/subtraction quadrature was not resolved.')
        mesh.append(dict(plasma_soft_cap=L,vacuum_and_thermal_residual=at.tolist(),
                         relative_to_frozen_born_before_gauge=float(at.sum()/stim),
                         comparison_quadrature_abs_defect=float(max(abs(at-small))),
                         cutoff_1e10_abs_defect=float(max(abs(at-mini))),
                         emission_angular_fraction=float(F(M,m,s,ep,q,L)[1]),
                         absorption_angular_fraction=float(F(M,m,s,ep,q,-L)[1])))
    return dict(parent_energy=float(ep),net_born_pauli=stim,zero_soft_radiator=float(y0),independent_angular_chart_defect=chart_error,
                angular_manufactured_soft_coefficients=c.tolist(),series_relative_defect=cd,
                first_kinematic_or_shell_knot=float(min(limits)),rows=mesh)


def run():
    rows=[]
    for grid in ('fine','third'):
        with gzip.open(base/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        _,lo,hi,*_=intervals(M,s,m,Q)
        rows.append(dict(grid=grid,data=[dict(window_fraction=f,result=one(d,M,m,s,lo+f*(hi-lo),Q))
                                        for f in (.2,.5,.8)]))
    return dict(scope='Three frozen Born parents at one scalar q, exact signed real recoil and leading virtual only, isotropic plasma photon bath.',rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
