"""Full frozen bare complex angular overlap; no physical damping inferred."""
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_complex_guard_reference import chart,local_cutoff_conversion
from finite_mass_qed_charged_closed_radial_pv_reference import angular_kernel,full_radial
from finite_mass_qed_charged_overlap_asymptotic_reference import solve_axes

D=np.longdouble
CD=np.clongdouble


def imag_conversion(v, cut):
    def H(h):
        h=CD(h)
        a,b=h+np.sqrt(1+h*h),h-np.sqrt(1+h*h)
        if abs(a)>abs(b): b=-1/a
        else:a=-1/b
        return -2*(np.log(-a)-np.log(-b))/(a-b)
    limit=H(D(v)-1j*D('1e-26')).imag
    def fn(log_t):
        return float(H(D(v)-1j*np.exp(-D(log_t))).imag/limit-(1 if log_t>0 else 0))
    return fsum(quad(fn,a,b,epsabs=2e-11,epsrel=2e-11,limit=200)[0] for a,b in ((-cut,0),(0,cut)))


def audit(d,E,cut,lo,hi,is_edge):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    name='parent' if cut=='parent' else 'fermion'
    call,_,data,C,Z,p,w,ma,A,B=chart(M,m,s,E,cut)
    raw=angular_kernel(M,m,s,E,cut)
    _,O,_=solve_axes(M,m,s,E,cut)
    S=np.asarray(d[f'endpoint_{name}_occupations'])
    b=np.hypot(d[f'endpoint_{name}_momentum_edges'],ma)-w
    if is_edge:
        j=np.argmin(abs(b))
        if abs(b[j])>1e-14:raise ArithmeticError('Declared archive step was not pinned.')
        b[j]=0.
    i,j=(np.searchsorted(b,0.,side=side)-1 for side in ('left','right'))
    Jr=O*p/2*(S[j]-S[i])
    Ji=np.pi*Z*(S[j]+S[i])/(2*Q*p*p*np.hypot(A,B))
    if Ji<=0:raise ArithmeticError('The occupied charged aperture must have a positive bare absorptive log.')
    R0=3e-4
    C_r=local_cutoff_conversion(B/A,40)
    C_i=imag_conversion(B/A,53)
    if abs(C_i-imag_conversion(B/A,44))>3e-9:
        raise ArithmeticError('Imaginary cutoff conversion tail did not settle.')

    def quadrature(N,eta):
        x,wgt=roots_legendre(N)
        yp=np.pi*(x+1)/2
        xp=np.sin(yp/2)**2
        floor=1e-23 if eta else 1e-11
        limits=[lo-E,hi-E]
        edge_knots=[t+h for t in limits for h in (0.,-1e-5,1e-5,-1e-7,1e-7,-1e-9,1e-9)]
        close=np.geomspace(floor,R0,22)
        if eta:
            spread=eta/(2*Q*p*np.hypot(A,B))
            close=sorted({*close,*[spread*r for r in (.03,.3,1.,3.,30.) if floor<spread*r<R0]})
        tail=np.geomspace(R0,max(abs(b[0]),b[-1]),49)
        breaks=sorted({*edge_knots,*[s*t for s in (-1,1) for t in (*close,*tail)]})
        parts_r=[];parts_i=[]
        for a1,b1,f in zip(b[:-1],b[1:],S,strict=True):
            if not f: continue
            for a,c in ((a1,min(-floor,b1)),(max(floor,a1),b1)):
                if c<=a: continue
                knots=sorted({a,c,*[v for v in breaks if a<v<c]})
                for l,r in zip(knots[:-1],knots[1:],strict=True):
                    pos=l+(r-l)*xp
                    fac=wgt*np.sin(yp)*(r-l)*np.pi*f/4
                    re=[]; im=[]
                    for v in pos:
                        mid,width,um,u0,g,pr,k=data(v)
                        if eta:
                            val=call(v,eta,original=True)
                            re.append(float(val.real*k/2));im.append(float(val.imag*k/2))
                        else:
                            re.append(float(raw(v)*k/2))
                            val=-np.pi*(C-2*Z/(pr*np.hypot(um,width)))/abs(g) if abs(u0)<1 else 0.
                            im.append(float(val*k/2))
                    parts_r.append(fac@re);parts_i.append(fac@im)
        return np.array([fsum(parts_r),fsum(parts_i)])

    rr=[]
    for eta in (1e-9,1e-10,1e-11,1e-12,1e-13,1e-14):
        a=quadrature(34,eta)
        b1=quadrature(54,eta)
        err=float(max(abs(a-b1)))
        width=eta/(2*Q*p*A)
        finite=b1-np.array([Jr*(np.log(R0/width)+C_r),Ji*(np.log(R0/width)+C_i)])
        rr.append(dict(angular_invariant_eta=eta,raw_real_and_imag=b1.tolist(),
                       hard_equivalent_finite_real_and_imag=finite.tolist(),quadrature_error=err))
        if err>2e-5:
            raise ArithmeticError(f'Whole support complex radial quadrature did not resolve: {err}.')
    hard=quadrature(42,0.)
    clipped=hard-np.array([Jr,Ji])*np.log(R0/1e-11)
    full=full_radial(d,E,cut,lo,hi,is_edge)['step']['value']
    # Compare real parts against an independent per-cell odd extraction.
    if abs(clipped[0]-full)>2e-5:
        raise ArithmeticError(f'Sharp full real replay did not reproduce per-cell extraction: {clipped[0]}, {full}.')
    moved=np.array(rr[-1]['hard_equivalent_finite_real_and_imag'])
    defects=abs(moved-clipped)
    slopes=(np.array(rr[-1]['raw_real_and_imag'])-rr[-2]['raw_real_and_imag'])/np.log(10.)
    if max(defects)>2e-4 or max(abs(slopes-[Jr,Ji]))>5e-5:
        raise ArithmeticError(f'Expected two cutoff logarithms did not close: {defects}, {slopes}, {Jr,Ji}.')
    return dict(external_parent_energy=E,cut=cut,tuned_to_archive_edge=is_edge,
                raw_log_coefficient_real_and_imag=[float(Jr),float(Ji)],
                measured_log_slope_real_and_imag=slopes.tolist(),conversion_numbers_real_and_imag=[C_r,C_i],
                exact_angular_hard_finite_parts=clipped.tolist(),independent_real_per_cell=full,
                conversion_defects=defects.tolist(),ladder=rr)


def run():
    out=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        k='daughter' if grid=='fine' else 'parent'
        name,ma,shift=('fermion',m,Es) if k=='daughter' else ('parent',M,0.)
        e=np.hypot(d[f'endpoint_{name}_momentum_edges'][1:-1],ma)+shift
        idx=np.flatnonzero((e>lo+4e-4)&(e<hi-4e-4))
        track=idx[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[idx]))]
        out.append(dict(grid=grid,rows=[audit(d,(lo+hi)/2,k,lo,hi,False),audit(d,float(e[track]),k,lo,hi,True)]))
    return dict(scope='Bare absorptive and dispersive overlap with a mathematical charged denominator guard; not a medium spectral function or physical damping.',results=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
