"""Complex angular guard for the bare positive charged overlap."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_closed_radial_pv_reference import angular_kernel
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes

D=np.longdouble
CD=np.clongdouble


def complex_polar(mid,width,pm):
    # u=mid+width*sinh y, v=exp(y). Roots have product -1.
    # Log differences use the same half-plane at either v endpoint.
    h=CD(pm/width)
    sq=np.sqrt(1+h*h)
    a,b=h+sq,h-sq
    if abs(a)>=abs(b):
        b=-1/a
    else:
        a=-1/b
    l,r=np.exp(np.arcsinh(np.array([-1-mid,1-mid],dtype=D)/width))
    if abs(a-b)<1e-22:
        raise ArithmeticError('Coalesced complex transformed roots.')
    val=(np.log(r-a)-np.log(l-a))-(np.log(r-b)-np.log(l-b))
    return 2*val/(width*(a-b))


def chart(M,m,s,E,species):
    M,m,s,E,q=map(D,(M,m,s,E,Q))
    es=np.hypot(s,q)
    z=(E*es-(M*M+s*s-m*m)/2)/q
    mass,w,a,sign=(M,E,z,D(1)) if species=='parent' else (m,E-es,z-q,D(-1))
    p=np.sqrt((w-mass)*(w+mass))
    c=a/p
    d=(M*M+m*m-s*s)/2
    C=(M+m)**2/(d+M*m)
    K=4-2*(M+m)*mass/(d+M*m)

    def data(delta):
        delta=D(delta)
        t=w+delta
        r=np.sqrt((t-mass)*(t+mass))
        R=p*r
        v=(t+w)/(r+p)
        dv=mass**2*(1/(t+r)+1/(w+p))/(r+p)
        b=delta**2*dv*(v+1)/2
        mid=c*(1+b/R)
        width=np.sqrt((1-c)*(1+c)*b*(2*R+b))/R
        um=delta*(es/q-c*v)/r-c*b/R
        u0=mid+um
        gp=2*sign*q*r
        return mid,width,um,u0,gp,R,r

    def evaluate(delta,eta,original=False):
        mid,width,um,u0,gp,R,r=data(delta)
        Z=CD(u0)-1j*D(eta)/gp
        zm=CD(um)-1j*D(eta)/gp
        v=complex_polar(mid,width,zm)
        l,r=np.arcsinh(np.array([-1-mid,1-mid],dtype=D)/width)
        ordinary=-(r-l)/(2*R)
        phase=np.log(CD(1)-Z)-np.log(CD(-1)-Z)
        val=(C*phase-2*d*v/R)/gp+K*ordinary/2
        if original:
            val+=K*1j*D(eta)*v/(4*R*gp)
        return val

    def independent(delta,eta,original):
        mid,width,um,u0,gp,R,r=map(float,data(delta))
        den=eta/gp
        def fun(u):
            G=gp*(u-u0)
            photon=-1/(2*R*np.hypot(u-mid,width))
            term=(C+4*d*photon)/(G+1j*eta)
            rem=K*photon/2
            if original:rem*=G/(G+1j*eta)
            return term+rem
        p=max(width,abs(u0-mid))
        guides=sorted({-1.,1.,mid,u0,*[mid+i*p for i in (-8,-1,1,8)],
                       *[u0+i*abs(den) for i in (-12,-3,-1,1,3,12)]})
        guides=[z for z in guides if -1<=z<=1]
        return (quad(lambda u:float(fun(u).real),-1,1,points=guides[1:-1],epsabs=1e-8,epsrel=1e-10,limit=330)[0]
                +1j*quad(lambda u:float(fun(u).imag),-1,1,points=guides[1:-1],epsabs=1e-8,epsrel=1e-10,limit=330)[0])
    A=mass*np.sqrt((1-c)*(1+c))/(p*p)
    B=(es/q-c*w/p)/p
    return evaluate,independent,data,float(C),float(d),float(p),float(w),float(mass),float(A),float(B)


def local_cutoff_conversion(b, resolution):
    """Offset from a hard radial clip to a smooth charged angle guard."""
    base=D(b)
    def H(h):
        h=CD(h)
        a,b=h+np.sqrt(1+h*h),h-np.sqrt(1+h*h)
        if abs(a)>abs(b):
            b=-1/a
        else:
            a=-1/b
        return -2*(np.log(-a)-np.log(-b))/(a-b)
    limit=H(base-1j*D("1e-25")).real
    def fn(log_t):
        t=np.exp(D(log_t))
        return float(H(base-1j/t).real/limit-(1 if log_t>0 else 0))
    return fsum(quad(fn,l,r,epsabs=2e-11,epsrel=1e-11,limit=150)[0]
                for l,r in ((-resolution,0),(0,resolution)))


def probe(d,E,species,edge):
    M=d['x_end'];m=d['mass_ratio']*M;s=d['checkpoints'][-1]['scalar_pole_over_T']
    K,ind,data,C,Ddot,p,w,mass,A,B=chart(M,m,s,E,species)
    bare=angular_kernel(M,m,s,E,species)
    _,odd,_=solve_axes(M,m,s,E,species)
    name='parent' if species=='parent' else 'fermion'
    H=np.hypot(d[f'endpoint_{name}_momentum_edges'],mass)-w
    occ=np.asarray(d[f'endpoint_{name}_occupations'])
    if edge:
        idx=np.argmin(abs(H))
        H[idx]=0.
        if not 0<idx<len(H)-1:raise ArithmeticError('Selected edge is outside archive cells.')
        jump=occ[idx]-occ[idx-1]
    else:
        jump=0.
    coeff=odd*p*jump/2
    nodes=[]
    for dr,width in ((-.00013,1e-8),(.00019,1.7e-8),(.00028,6.0e-8)):
        for scheme in (False,True):
            z=K(dr,width,scheme)
            alt=ind(dr,width,scheme)
            nodes.append(abs(z-alt)/max(1.,abs(z),abs(alt)))
    # The real part tends to the Cauchy PV; imaginary part to the delta
    # contribution with sign fixed by +i*eta in the charged denominator.
    lim=[]
    for dr in (-1e-4,9e-5):
        mid,wi,mu,u0,gp,R,_=data(dr)
        real=bare(dr)
        imag=-np.pi*(C-2*Ddot/(R*np.hypot(mu,wi)))/abs(gp)
        value=K(dr,1e-17)
        lim.append(max(abs(value.real-real)/max(1.,abs(real)),
                       abs(value.imag-imag)/max(1.,abs(imag))))
    if max(nodes)>2e-8 or max(lim)>5e-7:
        raise ArithmeticError(f'Complex angular branch or zero-width delta failed: {nodes}, {lim}.')
    def radial(n,eta):
        t,wt=roots_legendre(n)
        L=3e-4
        # The smoothed kernel is integrable at the exact radial collision.
        # Resolve the width in delta; omitting only delta<1e-14 makes a
        # negligible contribution in the tested guard hierarchy.
        stops=sorted({1e-14,L,*[x for x in abs(H) if 1e-14<x<L],
                      *np.geomspace(1e-14,L,19)})
        def val(x):
            rows=[]
            for scheme in (False,True):
                plus=H.searchsorted(x,side='right')-1
                minus=H.searchsorted(-x,side='right')-1
                left=K(x,eta,scheme).real if eta else bare(x)
                right=K(-x,eta,scheme).real if eta else bare(-x)
                a=left*data(x)[-1]*occ[plus]/2
                b=right*data(-x)[-1]*occ[minus]/2
                rows.append(float((a+b)*x))
            return rows
        parts=[]
        for l,r in zip(stops[:-1],stops[1:],strict=True):
            a,b=np.log(l),np.log(r)
            x=np.exp((a+b+(b-a)*t)/2)
            parts.append((b-a)/2*(wt@np.asarray([val(z) for z in x])))
        return np.sum(parts,axis=0)
    rows=[]
    for eta in (1e-9,1e-10,1e-11,1e-12,1e-13,1e-14):
        a=radial(12,eta)
        b=radial(20,eta)
        mesh=max(abs(a-b))
        if mesh>2e-6:
            raise ArithmeticError(f'Complex-guard radial mesh did not resolve: {a}, {b}.')
        rows.append(dict(imaginary_fermion_invariant=eta,
                         partial_fraction_pole_guard=b[0],original_charged_denominator_guard=b[1],
                         mesh_defect=float(mesh)))
    measured=(rows[-1]['partial_fraction_pole_guard']-rows[-2]['partial_fraction_pole_guard'])/np.log(10)
    limit=abs(measured-coeff)
    if limit>2e-4*max(1.,abs(coeff)):
        raise ArithmeticError(f'Smooth complex regulator failed its overlap log slope: {measured}, {coeff}.')
    K0=local_cutoff_conversion(B/A,32.)
    K1=local_cutoff_conversion(B/A,40.)
    if abs(K0-K1)>2e-9:
        raise ArithmeticError('Independent extent did not resolve the local log conversion.')
    L=3e-4
    hard_guard=1e-14
    hard=radial(20,0.)[0]-coeff*np.log(L/hard_guard)
    eta=rows[-1]['imaginary_fermion_invariant']
    delta0=eta/(2*Q*p*A)
    smooth=rows[-1]['partial_fraction_pole_guard']-coeff*(np.log(L/delta0)+K1)
    mismatch=abs(hard-smooth)
    if mismatch>5e-5:
        raise ArithmeticError(f'Analytic local conversion did not match the two radial prescriptions: {hard}, {smooth}, {K1}.')
    return dict(species=species,energy=E,tuned_to_archive_edge=edge,occupation_jump=float(jump),
                predicted_radial_log_coefficient=float(coeff),computed_cutoff_log_slope=float(measured),
                independent_angular_defect=float(max(nodes)),zero_width_pv_and_delta_defect=float(max(lim)),
                regulator_ladder=rows,log_conversion_from_local_kernel=K1,
                radial_width_for_last_smooth_guard=float(delta0),
                hard_subtracted_finite_part=float(hard),converted_smooth_finite_part=float(smooth),
                conversion_defect=float(mismatch))


def run():
    out=[]
    for n in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{n}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        k='daughter' if n=='fine' else 'parent'
        name,ma,add=('fermion',m,Es) if k=='daughter' else ('parent',M,0.)
        edges=np.hypot(d[f'endpoint_{name}_momentum_edges'][1:-1],ma)+add
        j=np.flatnonzero((edges>lo+4e-4)&(edges<hi-4e-4))
        ix=j[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[j]))]
        out.append(dict(grid=n,controls=[probe(d,(lo+hi)/2,k,False),probe(d,float(edges[ix]),k,True)]))
    return dict(scope='Mathematical complex guard on the bare positive charged angular pole; not a causal width, plasmon kernel or physical medium response.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
