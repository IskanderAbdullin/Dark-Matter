"""Charge-symmetric on-shell dispersive bridge, not a dressed decay rate."""
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import expit,roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals
from finite_mass_qed_bare_external_self_energy_reference import charged
from finite_mass_qed_log_thermodynamic_reference import shell_integral

D=np.longdouble


def bare_branches(p,m,k):
    # Independent signed bare angular kernel, evaluated in k rather
    # than by subtracting two nearly equal energies.
    P,M,K=map(D,(p,m,k))
    E,e=np.hypot(P,M),np.hypot(K,M)
    if not P:
        g=np.array([2*M*(M+e),-2*M*K*K/(M+e)],dtype=D)
        return (2+4*M*M/g)/E
    R=P*K
    diff=(K-P)*(K+P)/(e+E)
    pos_hi=E*e-M*M+R
    pos_lo=M*M*diff*diff/pos_hi
    neg_hi=E*e+M*M+R
    neg_lo=M*M*(E+e)**2/neg_hi
    if R:
        jp=np.log(pos_lo/pos_hi)/(2*R) if pos_lo<pos_hi*.3 else np.log1p(-2*R/pos_hi)/(2*R)
        jn=np.log(neg_hi/neg_lo)/(2*R) if neg_lo<neg_hi*.3 else np.log1p(2*R/neg_lo)/(2*R)
    else:
        jp,jn=-1/pos_hi,1/neg_hi
    return (2+2*M*M*np.array([jn,jp]))/E


def cuts(p,m,edges,weights,occ_fn=None,N=36):
    nodes,nw=roots_legendre(N)
    vals=[]
    for lo,hi,f in zip(edges[:-1],edges[1:],weights,strict=True):
        if not f:
            continue
        split=sorted({lo,hi,*[x for x in (p,p-.0002,p+.0002,p-.02,p+.02, m) if lo<x<hi]})
        for a,b in zip(split[:-1],split[1:],strict=True):
            z=np.pi*(nodes+1)/2
            ks=a+(b-a)*np.sin(z/2)**2
            jac=(b-a)*np.pi/4*np.sin(z)
            e=np.hypot(ks,m)
            w=ks*ks/e/2*jac*nw*f
            if occ_fn is not None:
                w*=occ_fn(ks)
            ys=np.array([bare_branches(p,m,k) for k in ks])
            vals.append(w@ys)
    return np.sum(vals,axis=0).astype(float)


def equilibrium(p,m):
    E=np.hypot(p,m)
    bounds=np.r_[0.,.008,.03,.12,.4,.9,1.8,4,8,16,26,46.]
    fd=lambda k:expit(-np.hypot(k,m))
    a=cuts(p,m,bounds,np.ones(len(bounds)-1),fd,29)
    b=cuts(p,m,bounds,np.ones(len(bounds)-1),fd,52)
    j=quad(lambda k: k*k*fd(k)/np.hypot(k,m),0,46,epsabs=7e-13,epsrel=2e-12)[0]
    log=shell_integral(p,m,46.)
    db=2*(np.pi**2/6)/E
    expected=(2*j-m*m*log)/E+db
    actual=fsum(b)+db
    if abs(actual-expected)>3e-9 or max(abs(a-b))>3e-8:
        raise ArithmeticError(f'Equilibrium pole and bare cut shapes do not agree: {p,m,a,b,actual,expected}.')
    return dict(m_over_T=m,p_over_T=p,charged_negative_and_positive=b.tolist(),
                photon_with_thermal_sign=db,bare_energy_shift_before_gauge=actual,
                published_equilibrium_formula_before_gauge=expected,
                absolute_defect=abs(actual-expected),numerical_bare_refinement=float(max(abs(a-b))))


def archived_log(p,m,edges,occ):
    E=np.hypot(p,m)
    fac=sing=0.
    for l,h,f in zip(edges[:-1],edges[1:],occ,strict=True):
        if not f:
            continue
        j=quad(lambda k:k*k/np.hypot(k,m),l,h,epsabs=2e-13)[0]
        fac+=f*j
        if p==0:
            r=2*quad(lambda k:1/np.hypot(k,m),l,h,epsabs=2e-13)[0]
        else:
            def A(k):
                return k/np.hypot(k,m)
            z=quad(lambda k:A(k)*np.log(p+k),l,h,epsabs=2e-12)[0]
            if l<=p<=h:
                for sign,limit in ((-1,p-l),(1,h-p)):
                    if limit>0:
                        z-=quad(lambda x:A(p+sign*x),0,limit,weight='alg-loga',wvar=(0.,0.),epsabs=3e-13)[0]
            else:
                z-=quad(lambda k:A(k)*np.log(abs(p-k)),l,h,epsabs=2e-12)[0]
            r=z/p
        sing+=f*r
    return np.array([2*fac/E,-m*m*sing/E])


def archive(d):
    M=d['x_end'];m=M*d['mass_ratio']
    ms=d['checkpoints'][-1]['scalar_pole_over_T']
    es,lo,hi,*_=intervals(M,ms,m,Q)
    E=(lo+hi)/2
    checks=[]
    for name,ma,ext in (('parent',M,E),('fermion',m,E-es)):
        ed=np.asarray(d[f'endpoint_{name}_momentum_edges'])
        oc=np.asarray(d[f'endpoint_{name}_occupations'])
        p0=np.sqrt((ext-ma)*(ext+ma))
        for p in (0.,.03,.25,p0,4.):
            a=cuts(p,ma,ed,oc,N=38)
            b=cuts(p,ma,ed,oc,N=58)
            c=archived_log(p,ma,ed,oc)
            defect=abs(fsum(b)-fsum(c))
            if defect>2e-8 or max(abs(a-b))>2e-8:
                raise ArithmeticError(f'Archived cut energy bridge was not resolved: {a,b,c}.')
            if p>0:
                replay=[]
                for k in (.017,.4,2.4):
                    r=np.hypot(ma,k)-np.hypot(ma,p)
                    C=np.array([charged(ma,np.hypot(ma,p),r,sg)[0] for sg in (-1,1)])
                    replay.append(max(abs(C-bare_branches(p,ma,k))))
                if max(replay)>2e-9:
                    raise ArithmeticError('Independently coded pointwise shell did not replay the published derivative kernel.')
            checks.append(dict(species=name,external_bare_p=p,event_momentum=bool(p==p0),
                               signed_charged_projections=b.tolist(),factorizable_and_log=c.tolist(),
                               relative_log_over_factorized=float(c[1]/c[0]),
                               absolute_two_chart_defect=defect))
    return checks


def run():
    eq=[equilibrium(p,m) for m in (.003,.05,.25,1.,3.) for p in (0.,.017,.2,1.4,4.)]
    a=[]
    for n in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{n}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f: d=json.load(f)
        a.append(dict(grid=n,checks=archive(d)))
    return dict(scope='Charge-symmetric free-line on-shell dispersive comparison with independently known equilibrium pole result. Archive checks do not compute physical retarded residues or a matched medium rate.',
                bare_units='multiply shifts by g²/(4 pi²); archive includes no photon occupation assumed',
                equilibrium=eq,archive=a)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
