"""Gain/loss compatibility gate for the shared positive charged-cut support."""
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

from early_pole_decay_mode_clock_reference import Q, intervals


def setup(d, q):
    M = d['x_end']
    m = M*d['mass_ratio']
    s = d['checkpoints'][-1]['scalar_pole_over_T']
    E, lo, hi, span, _ = intervals(M, s, m, q)
    A = (M*M-m*m+s*s)/2
    Z = (M*M+m*m-s*s)/2
    H = np.sqrt((A-M*s)*(A+M*s))
    c = -np.pi*(M+m)**2/(4*q*(Z+M*m))
    k = np.pi*Z/(2*H)
    ep = np.hypot(M, d['endpoint_parent_momentum_edges'])
    eb = np.hypot(m, d['endpoint_fermion_momentum_edges'])+E
    fp = np.asarray(d['endpoint_parent_occupations'])
    fb = np.asarray(d['endpoint_fermion_occupations'])
    return E, lo, hi, span, c, k, ep, eb, fp, fb


def products(p, b):
    return p*(1-b), b*(1-p)


def thermal(E, lo, hi, c, k, center):
    # Manufactured common-temperature FD legs; only the shared statistical
    # factors are under test, with no radiative coupling or contour assigned.
    out = []
    for T, mu_b, x_goal in ((1., .12, 0.7), (.04, -.015, 1.4),
                             (.3, .4, 2e-5), (1., 0., E)):
        mu_p = mu_b+E-T*x_goal
        x = (E-(mu_p-mu_b))/T
        bose = 1/np.expm1(x)
        points = np.linspace(lo, hi, 201)
        p, b = expit(-(points-mu_p)/T), expit(-(points-E-mu_b)/T)
        G, L = products(p, b)
        rel = float(np.max(abs(G-np.exp(-x)*L)/L))
        net = float(np.max(abs((1+bose)*G-bose*L)/(bose*L)))
        # Two independent positive-cut radial weights at fixed guard:
        # the uniform Born measure and the full common absorptive kernel.
        def weight_integral(fn, eps):
            total = []
            for side, end in ((-1., center-lo), (1., hi-center)):
                a, b = np.log(eps), np.log(end)
                total.append(quad(lambda u: (c*np.exp(u)+k)*fn(center+side*np.exp(u)),
                                  a,b,epsabs=1e-13,epsrel=3e-11)[0])
            return fsum(total)
        def at(y, i):
            return products(expit(-(y-mu_p)/T),expit(-(y-E-mu_b)/T))[i]
        guards = []
        for eps in (1e-6, 1e-9, 1e-12):
            g, l = (weight_integral(lambda y:at(y,i),eps) for i in (0,1))
            b_g,b_l = (quad(lambda y:at(y,i),lo,hi,epsabs=1e-13,epsrel=3e-11)[0]
                       for i in (0,1))
            kern = float(np.log1p((l-g)/g))
            born = float(np.log1p((b_l-b_g)/b_g))
            signed = float(abs(-(l-g)-(-np.expm1(x)*g))/l)
            agree = max(abs(kern-x), abs(born-x), signed)
            if agree > 8e-14:
                raise ArithmeticError(('common-temperature guard ratio', eps, x, agree))
            guards.append(dict(guard=eps, log_loss_gain_signed_kernel=kern,
                               log_loss_gain_Born=born, signed_kernel_gain=g, signed_kernel_loss=l,
                               signed_difference_identity_error=signed))
        if max(rel,net)>8e-14:
            raise ArithmeticError(('pointwise KMS-style balance fixture', rel, net))
        out.append(dict(T=T, mu_parent=mu_p, mu_daughter=mu_b, born_affinity=x,
                        equilibrium_bose=bose, max_relative_pointwise_product_error=rel,
                        max_relative_pointwise_gain_loss_error=net, guards=guards))
    return out


def step_archive(lo,hi,c,k,ep,eb,fp,fb,center,tuned=None):
    gp,gb=ep-center,eb-center
    if tuned is not None:
        a={'parent':gp,'daughter':gb}[tuned]
        j=np.argmin(abs(a))
        if abs(a[j])>1e-13:
            raise ArithmeticError(('unresolved tuned edge',a[j]))
        a[j]=0.
    bounds = lo-center,hi-center
    turn = -k/c
    knots = sorted({*bounds,0.,*[p for p in (*gp,*gb,-turn,turn) if bounds[0]<p<bounds[1]]})
    def read(g,f,z):
        j=np.searchsorted(g,z,side='right')-1
        return f[j] if 0<=j<len(f) else 0.
    def parts(x):
        p,b=read(gp,fp,x),read(gb,fb,x)
        return *products(p,b),p-b
    cells=[]
    for a,b in zip(knots[:-1],knots[1:],strict=True):
        g,l,s=parts((a+b)/2)
        if abs((g-l)-s)>8e-17:
            raise ArithmeticError('archive Pauli composition did not close')
        cells.append((a,b,g,l,s))

    def weight(a,b,eps,kind):
        if eps is not None:
            a,b=(a,min(b,-eps)) if b<=0. else (max(a,eps),b)
        if b<=a:
            return 0.
        if kind=='uniform':
            return b-a
        small,large=sorted((abs(a),abs(b)))
        # Accurately integrate the kernel on narrow cells.
        return c*(b-a)+k*np.log1p((large-small)/small)

    def rates(eps,kind):
        items=[(weight(a,b,eps,kind),g,l,s) for a,b,g,l,s in cells]
        J=[fsum(i[0]*i[j] for i in items) for j in (1,2,3)]
        if abs((J[0]-J[1])-J[2])>3e-13*fsum(abs(i[0])*(i[1]+i[2]) for i in items):
            raise ArithmeticError('Weighted archive Pauli identity did not close')
        return J

    # Pointwise neutrality requires the same scalar occupation in every
    # cell. A moment-balanced weight can hide compensating residuals.
    def diagnose(r,z,kind):
        g,l,s=r
        bose=g/(l-g) if l>g and g>0 else None
        x=float(np.log(l/g)) if l>0 and g>0 else None
        if bose is None:
            bound=(g-l)/(g+l) if g>0 and 0<=l<g and kind=='uniform' else None
            if bound is not None:
                for n in (0.,1.,100.,1e8):
                    ratio=(g+n*(g-l))/(g+n*(g+l))
                    if not bound-1e-15<=ratio<=1.+1e-15:
                        raise ArithmeticError('Positive Born measure violated the neutrality bound')
            return dict(weight=kind,guard=z,gain=g,loss=l,signed_difference=s,
                        affinity=x,neutral_bose=None,
                        all_n_ge_zero_net_over_turnover_lower_bound=bound)
        terms=[((1+bose)*G-bose*L,weight(a,b,z,kind),(1+bose)*G+bose*L)
               for a,b,G,L,_ in cells]
        cancellation=fsum(u*w for u,w,_ in terms)
        turnover=fsum(abs(w)*t for _,w,t in terms)
        if abs(cancellation)>5e-13*turnover:
            raise ArithmeticError('Balanced mean was not neutral')
        return dict(weight=kind,guard=z,gain=g,loss=l,signed_difference=s,
                    affinity=x,neutral_bose=bose, signed_net_over_turnover=cancellation/turnover,
                    gross_local_net_over_turnover=fsum(abs(u*w) for u,w,_ in terms)/turnover,
                    negative_weight_fraction_of_abs_turnover=fsum(-w*t for _,w,t in terms if w<0)/turnover,
                    number_of_positive_cells=int(sum(u*w>0 for u,w,_ in terms)),
                    number_of_negative_cells=sum(int(u*w<0) for u,w,_ in terms))
    report=[diagnose(rates(None,'uniform'),None,'uniform')]
    for eps in (1e-5,1e-7,1e-9,1e-12):
        report.append(diagnose(rates(eps,'common'),eps,'common'))
    left=next(z for z in cells if z[1]==0.)
    right=next(z for z in cells if z[0]==0.)
    G=left[2]+right[2]
    L=left[3]+right[3]
    alpha=float(np.log(L/G)) if L>0 and G>0 else None
    # Exact expansion of all three signed/unsigned profiles after guard
    # falls into the two central cells. Their constants cancel in
    # differences between the two guards.
    e1=min(1e-9,-left[0]/2,right[1]/2)
    e2=e1/1000
    R1=np.array(rates(e1,'common'))
    R2=np.array(rates(e2,'common'))
    pred=np.array(left[2:])+np.array(right[2:])
    error=float(np.max(abs(R2-R1-pred*(k*np.log(e1/e2)+c*(e1-e2)))))
    if error>8e-13:
        raise ArithmeticError(('two-trace radial guard expansion',error))
    # Zero crossings extrapolated analytically within fixed step cells.
    # Deliberately return log10 guard, which can be far below floating-point.
    crossing=[]
    for v,slope in zip(R1[:2],pred[:2],strict=True):
        t=-v/(k*slope)-c*e1/k if v<0 and slope>0 else None
        crossing.append(float(np.log10(e1)-t/np.log(10.)) if t is not None and t>50 else None)

    # Direct log-radial quadrature, without the integrated primitive.
    independent=[]
    for a,b,g,l,sgn in cells:
        start,finish=sorted((abs(a),abs(b)))
        start=max(start,1e-9)
        if finish>start:
            w=quad(lambda u:c*np.exp(u)+k,np.log(start),np.log(finish),
                   epsabs=1e-9,epsrel=2e-12)[0]
            independent.append(w*np.array([g,l,sgn]))
    replay=float(np.max(abs(np.sum(independent,axis=0)-R1))) if e1==1e-9 else float(
        np.max(abs(np.sum(independent,axis=0)-np.array(rates(1e-9,'common')))))
    if replay>1e-9:
        raise ArithmeticError(('independent stepped archive log quadrature',replay))
    return dict(external_parent_E=center,tuned=tuned,cells=len(cells),
                kernel_sign_flip_radial_distance=turn,
                guard_terms_trace_G_left_right=[left[2],right[2]],
                guard_terms_trace_L_left_right=[left[3],right[3]],
                zero_guard_log_loss_gain_limit=alpha,guard_difference_abs_defect=error,
                extrapolated_log10_guard_raw_gain_and_loss_zeros=crossing,
                independent_signed_archive_quadrature_abs_defect=replay,
                log_loss_gain_supported_cell_range=[min(np.log(L/G) for _,_,G,L,_ in cells if L>0 and G>0),
                                                    max(np.log(L/G) for _,_,G,L,_ in cells if L>0 and G>0)],
                gain_zero_only_cells=sum(int(g==0 and l>0) for _,_,g,l,_ in cells),
                loss_zero_only_cells=sum(int(l==0 and g>0) for _,_,g,l,_ in cells),weights=report)


def run():
    data=[]
    for grid in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        E,lo,hi,span,c,k,ep,eb,fp,fb=setup(d,Q)
        mids=[lo+f*span for f in (.20,.5,.8)]
        rows=[step_archive(lo,hi,c,k,ep,eb,fp,fb,u) for u in mids]
        for tag,edge,occ in [('parent',ep,fp),('daughter',eb,fb)]:
            candidates=np.flatnonzero((edge[1:-1]>lo+4e-4)&(edge[1:-1]<hi-4e-4))
            if len(candidates):
                i=candidates[np.argmax(abs(np.diff(occ)[candidates]))]
                rows.append(step_archive(lo,hi,c,k,ep,eb,fp,fb,edge[i+1],tag))
        data.append(dict(grid=grid,q=Q,scalar_energy=E,parent_edges=[lo,hi],
                         common_constant=c,common_inverse_coefficient=k,
                         thermal=thermal(E,lo,hi,c,k,mids[1]),archive=rows))
    return dict(scope='Bare shared positive charged-cut support and kernel: algebraic gain/loss compatibility only; not an assignment of thermal vertex, cut rate, or physical scalar occupation.',data=data)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
