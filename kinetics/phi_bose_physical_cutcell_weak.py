"""Independent physical Yukawa weak source integrated over exact daughter cut cells.

The first algorithm splits parent momentum at *on-shell cell crossings* and
integrates Fermi/Bose gain-loss within each cell with positive log-energy Gauss
weights. A second algorithm integrates exact FD*BE daughter log/dilog primitives
by independent adaptive parent quadrature. The former is NOT an evolving
nonlinear strong collision operator or a relic prediction.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, spence

M1, M2, MPHI, YUKAWA, G1 = .1, .02, 86e-9, .001, 4
LAM = (M1*M1-(M2+MPHI)**2)*(M1*M1-(M2-MPHI)**2)
KSTAR = math.sqrt(LAM)/(2*M1)
ESTAR = math.hypot(KSTAR,MPHI)
WIDTH = YUKAWA**2*((M1+M2)**2-MPHI**2)*math.sqrt(LAM)/(16*math.pi*M1**3)


def shell(p):
    e1=math.hypot(p,M1)
    b=p*KSTAR/M1
    kmin=(p*p*MPHI*MPHI+M1*M1*ESTAR*ESTAR)/(M1*(e1*ESTAR+p*KSTAR))
    kmax=(e1*ESTAR+p*KSTAR)/M1
    return e1,kmin,kmax,b


def scalar_edges(pmax,n):
    upper=math.hypot(pmax,M1)+.2
    return MPHI+(upper-MPHI)*(np.arange(n+1,dtype=float)/n)**2


def scalar_capacity(lo,hi):
    """Physical single-state number capacity V and positive radial energy W."""
    pl=math.sqrt(max(0,lo*lo-MPHI*MPHI))
    ph=math.sqrt(max(0,hi*hi-MPHI*MPHI))
    v=(ph**3-pl**3)/(6*math.pi**2)
    x,w=leggauss(48)
    ps=pl+(ph-pl)*(x+1)/2
    q=(ph-pl)*w/2
    energy=math.fsum(q*ps*ps*np.sqrt(ps*ps+MPHI*MPHI))/(2*math.pi**2)
    return v,energy


def parent_splits(pmax,edges):
    """Partition at physical k_min/k_max intersections with each cell edge."""
    cuts=[0.,pmax]
    kmin_pmax=shell(pmax)[1]
    kmax_pmax=shell(pmax)[2]
    for edge in edges[1:-1]:
        if kmin_pmax < edge < ESTAR:
            cuts.append(brentq(lambda p: shell(p)[1]-edge,0.,pmax,xtol=1e-12))
        elif ESTAR < edge < kmax_pmax:
            cuts.append(brentq(lambda p: shell(p)[2]-edge,0.,pmax,xtol=1e-12))
    return sorted(cuts)


def cutcell_weak(pmax,nrad,nangle,ncell,t1=.9,t0=1.,mu1=0.,mu2=0.,
                 muphi=0.,split_parent=True):
    """Return Rphi, Qphi, Q1, Q2 and physically positive cell capacities.

    nrad is Gauss order PER smooth parent segment if split_parent=True.
    """
    assert pmax>0 and nrad>0 and nangle>0 and ncell>=2
    assert t1>0 and t0>0 and muphi<MPHI
    edges=scalar_edges(pmax,ncell)
    x,w=leggauss(nrad)
    z,v=leggauss(nangle)
    segments=parent_splits(pmax,edges) if split_parent else [0.,pmax]
    bycell=np.zeros((ncell,4),float)  # Rphi, Qphi, Q1, Q2
    heat=0.; minpos=math.inf
    for left,right in zip(segments[:-1],segments[1:]):
        for px,pw in zip(x,w):
            p=left+(right-left)*(px+1)/2
            dp=(right-left)*pw/2
            e1,kmin,kmax,b=shell(p)
            parent_weight=G1*p*p/(2*math.pi**2)*WIDTH*M1/e1*dp
            f1=expit(-(e1-mu1)/t1)
            for j,(low,high) in enumerate(zip(edges[:-1],edges[1:])):
                l=max(low,kmin);h=min(high,kmax)
                if h<=l: continue
                d=math.log(h/l)
                k=l*np.exp((z+1)*d/2)
                dk=k*d*v/2
                e2=e1-k
                f2=expit(-(e2-mu2)/t0)
                bose=1/np.expm1((k-muphi)/t0)
                forward=f1*(1-f2)*(1+bose)
                reverse=(1-f1)*f2*bose
                net=forward-reverse
                W=parent_weight*dk/(2*b)
                flux=W*net
                n=float(np.sum(flux))
                qphi=float(flux@k)
                q2=float(flux@e2)
                bycell[j]+=(n,qphi,-e1*n,q2)
                affinity=-(e1-mu1)/t1+(e1-mu2-muphi)/t0
                heat+=n*affinity
                minpos=min(minpos,float(np.min(W)))
    total=np.sum(bycell,axis=0)
    capacities=np.array([scalar_capacity(l,h) for l,h
                         in zip(edges[:-1],edges[1:])])
    return dict(Rphi=float(total[0]),Qphi=float(total[1]),Q1=float(total[2]),
                Q2=float(total[3]),energy_defect=float(sum(total[1:])),
                entropy_rate=float(heat),bycell=bycell,capacities=capacities,
                edges=edges,parent_segments=len(segments)-1,
                min_event_weight=minpos)


def primitives(e1,lo,hi,t0,mu2,muphi):
    """Exact independent daughter FD*BE number and phi-energy primitives."""
    a=math.exp((e1-mu2-muphi)/t0)
    loga=math.log(a)
    def primitive(k):
        u=(k-muphi)/t0
        t=math.exp(u)
        L=math.log(-math.expm1(-u))+u-float(np.logaddexp(u,loga))
        P=-spence(t)-u*math.log1p(t/a)-spence(1+t/a)
        return L,P
    l0,l1=primitive(lo);h0,h1=primitive(hi)
    n=t0*(h0-l0)/(1+a)
    q=muphi*n+t0*t0*(h1-l1)/(1+a)
    return n,q


def analytic_weak(pmax,t1=.9,t0=1.,mu1=0.,mu2=0.,muphi=0.):
    """Independent adaptive 1D parent integral after exact daughter integration."""
    def f(p):
        e,lo,hi,b=shell(p)
        if b==0:return (0.,0.,0.)
        n,q=primitives(e,lo,hi,t0,mu2,muphi)
        f1=expit(-(e-mu1)/t1)
        aff=-(e-mu1)/t1+(e-mu2-muphi)/t0
        pref=G1*p*p/(2*math.pi**2)*WIDTH*M1/e*(1-f1)*math.expm1(aff)/(2*b)
        return pref*n,pref*q,-pref*e*n
    a=[quad(lambda p:f(p)[i],0.,pmax,epsabs=1e-19,
            epsrel=2e-10,limit=300)[0] for i in range(3)]
    return dict(Rphi=a[0],Qphi=a[1],Q1=a[2],Q2=-a[1]-a[2])


def audit():
    cases={'temperature_split':(.9,1.,0.,0.,0.),
           'chemical_forward':(1.,1.,-.09,.05,-.2),
           'chemical_reverse':(1.,1.,-.21,.05,-.2),
           'chemical_equilibrium':(1.,1.,-.15,.05,-.2)}
    rows=[]
    for case,settings in cases.items():
        for pmax in (2.,8.,16.):
            reference=analytic_weak(pmax,*settings)
            for nr,na,nc in ((6,8,8),(10,16,16),(16,32,32)):
                got=cutcell_weak(pmax,nr,na,nc,*settings)
                errors={k:abs(got[k]/reference[k]-1)
                        if abs(reference[k])>1e-22 else abs(got[k])
                        for k in ('Rphi','Qphi','Q1','Q2')}
                rows.append(dict(case=case,pmax=pmax,mesh=(nr,na,nc),
                                 reference=reference,relative_or_absolute_errors=errors,
                                 Rphi=got['Rphi'],Qphi=got['Qphi'],
                                 Q1=got['Q1'],Q2=got['Q2'],
                                 H=got['entropy_rate'],
                                 defect=got['energy_defect'],
                                 firstcell_R=got['bycell'][0,0],
                                 firstcell_Q=got['bycell'][0,1],
                                 segments=got['parent_segments']))
    return {'schema_version':1,'scope':'prescribed thermal weak Yukawa sources only, not strong kinetic evolution',
            'parameters_GeV':dict(M1=M1,M2=M2,MPHI=MPHI,YUKAWA=YUKAWA,G1=G1),
            'cases':rows,'class_ready':False}


if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2,allow_nan=False))
