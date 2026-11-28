"""Compare two mathematical regulators along a correlated overlap path."""
from functools import lru_cache
import gzip,json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q,intervals
from finite_mass_qed_charged_complex_guard_reference import chart

F=np.longdouble


@lru_cache(None)
def GL(n):
    return roots_legendre(n)


def polars(a,b,pole,core,width,n):
    """Geometrically resolve the photon peak and the displaced charged pole."""
    x,w=GL(n)
    radius=core*np.geomspace(1/256,1e36,118)
    charge=np.array([-24.,-7.,-2.,-.5,.5,2.,7.,24.])*width+pole
    limits=sorted({a,b,0.,pole,*charge,*radius,*(-radius)})
    limits=np.array([c for c in limits if a<=c<=b])
    l,h=limits[:-1],limits[1:]
    g=(h-l)/2
    q=l[:,None]+g[:,None]*(x+1)
    return q,g[:,None]*w


def audit(d,E,cut,step):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    _,_,dat,C,z,p,e,ma,A,B=chart(M,m,s,E,cut)
    es=np.hypot(F(s),F(Q))
    u=(F(E)*es-(F(M)**2+F(s)**2-F(m)**2)/2)/F(Q)
    if cut=='daughter':u-=Q
    c=float(u/p)
    sign=1 if cut=='parent' else -1
    R0=p*p
    den0=sign*2*Q*p
    axis=np.hypot(A,B)
    a,b=A/axis,B/axis
    name='parent' if cut=='parent' else 'fermion'
    span=np.hypot(d[f'endpoint_{name}_momentum_edges'],ma)-e
    occ=np.asarray(d[f'endpoint_{name}_occupations'])
    if step:
        ix=np.argmin(abs(span))
        if abs(span[ix])>5e-15:raise ArithmeticError('Declared step not on the archive bin.')
        span[ix]=0.
    side=[float(occ[np.searchsorted(span,0.,side=j)-1]) for j in ('left','right')]
    L=min(.0003,.55*min(abs(v) for v in span if abs(v)>1e-14))
    reg_piece=4-2*(M+m)*ma/(z+M*m)

    def integrate(eta,lam,n,mn,*,near=False,tail=None):
        radius=eta/(abs(den0)*axis)
        shift=axis*radius
        photon_width_unsigned=lam*R0*shift**2/(1-c*c)
        xmax=L/radius if tail is None else tail
        xmin=1e-9
        t,w=GL(mn)
        breaks=np.geomspace(xmin,xmax,32)
        total=np.zeros(2,dtype=complex)
        for low,high in zip(breaks[:-1],breaks[1:],strict=True):
            xl,xh=np.log(low),np.log(high)
            xs=np.exp(xl+(xh-xl)/2*(t+1))
            W=(xh-xl)/2*w
            for j,sg in enumerate((-1,1)):
                values=[]
                for x in xs:
                    xx=x*sg
                    phsign=sg*sign
                    if near:
                        photon=(a*xx)**2
                        dev=-1j*phsign*lam
                        yp=b*xx
                        y,V=polars(-1e13,1e13,yp,max(1e-40,min(abs(a*xx),np.sqrt(lam))),1.,n)
                        G0=np.sqrt(y*y+photon)
                        G1=np.sqrt(y*y+photon+dev)
                        phdiff=dev/(2*G0*G1*(G0+G1))
                        expr=4*z/(R0*den0*axis)*phdiff/(y-yp+1j)
                    else:
                        delta=xx*radius
                        mid,ww,um,u0,g,R,r=dat(delta)
                        R,g=float(R),float(g)
                        s0=float(ww)/shift
                        yp=float(um)/shift
                        ch=abs(eta/(g*shift))
                        y,V=polars(float((-1-mid)/shift),float((1-mid)/shift),yp,max(1e-40,min(s0,np.sqrt(lam))),ch,n)
                        realP=y*y+s0*s0
                        u=float(mid)+shift*y
                        # b_gap/R is quadratic in radial displacement. It
                        # matters only in the photon continuation correction.
                        dt=e+delta
                        rs=float(r)
                        v=(dt+e)/(rs+p)
                        b_gap=delta*delta*ma*ma*(1/(dt+rs)+1/(e+p))/(rs+p)*(v+1)/2
                        wp=phsign*photon_width_unsigned
                        dev=-1j*wp*(b_gap+R*(1-c*u))/(R*shift)**2-(wp/(2*R*shift))**2
                        G0=np.sqrt(realP)
                        G1=np.sqrt(realP+dev)
                        phdiff=dev/(2*G0*G1*(G0+G1))
                        ph=-1/(2*G1)
                        g0=g*shift*(y-yp)
                        expr=((4*z+reg_piece*g0/2)*phdiff-1j*C*wp*ph)/(R*g*axis)/(y-yp+1j*sign*eta/(g*shift))
                    values.append(np.sum(V*expr)*(p if near else float(r))/2*x)
                total[j]+=W@np.asarray(values)*side[j]
        return total

    rows=[]
    for lam in (.03,1.,30.):
        ref=integrate(1e-16,lam,18,10,near=True,tail=1e10)
        edged=integrate(1e-16,lam,18,10,near=True,tail=1e9)
        test=integrate(1e-16,lam,26,15,near=True,tail=1e10)
        if max(abs(ref-test))>8e-7 or max(abs(edged-test))>8e-7:
            raise ArithmeticError(f'Local two-width convergence failed: {cut}, {lam}, {ref}, {edged}, {test}.')
        sym=max(abs(test[0].real/side[0]+test[1].real/side[1]),
                abs(test[0].imag/side[0]-test[1].imag/side[1]))
        if sym>1e-8:
            raise ArithmeticError('Local two-width comparison lost odd real or even imaginary behavior.')
        vals=[]
        for eta in (1e-12,1e-14,1e-16):
            v=integrate(eta,lam,18,10)
            local=integrate(eta,lam,18,10,near=True)
            vals.append(dict(charged_eta=eta,exact_shift=[v.sum().real,v.sum().imag],
                             local_same_support=[local.sum().real,local.sum().imag],
                             difference_at_same_support=float(abs((v-local).sum()))))
        if vals[-1]['difference_at_same_support']>6e-5:
            raise ArithmeticError('Exact correlated regulator shift did not converge to its local model.')
        rows.append(dict(dimensionless_photon_width=lam,
                         local_asymptotic_shift_real_and_imag=[test.sum().real,test.sum().imag],
                         local_side_real=[z.real for z in test],local_side_imaginary=[z.imag for z in test],
                         parity_defect=float(sym),local_mesh_defect=float(max(abs(ref-test))),
                         local_tail_extent_defect=float(max(abs(edged-test))),ladder=vals))
    return dict(cut=cut,external_parent_energy=E,tuned_step=step,occupations_left_right=side,
                local_relative_angular_scales=[a,b],outer_half_width=L,comparators=rows)


def run():
    out=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        cut='daughter' if grid=='fine' else 'parent'
        name,ma,st=('fermion',m,Es) if cut=='daughter' else ('parent',M,0.)
        ee=np.hypot(d[f'endpoint_{name}_momentum_edges'][1:-1],ma)+st
        who=np.flatnonzero((ee>lo+4e-4)&(ee<hi-4e-4))
        j=who[np.argmax(abs(np.diff(d[f'endpoint_{name}_occupations'])[who]))]
        out.append(dict(grid=grid,rows=[audit(d,(lo+hi)/2,cut,False),audit(d,float(ee[j]),cut,True)]))
    return dict(scope='Regulator-path test for two bare scalar-latest charged overlap denominators, not a causal physical width or medium rate.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
