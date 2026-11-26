"""Discontinuous cell-wise entropy-dual Galerkin for on-shell Yukawa 1<->2.

Each physical energy cell has TWO independent endpoint inverse-fugacity
coordinates. Continuous within each cell, discontinuous across boundaries.
This preserves affine E and partition unity per collision but eliminates
cross-cell inverse susceptibility. This is a numerical IR realisability
pilot, NOT proof of dark matter condensation, abundance or FRW physics.
"""
from dataclasses import dataclass
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from scipy.linalg import solve


@dataclass(frozen=True)
class CellMode:
    mass: float
    edges: np.ndarray
    fermion: bool = False
    degeneracy: float = 1.

    def __post_init__(self):
        v = np.asarray(self.edges, dtype=float)
        if not (math.isfinite(self.mass) and self.mass >= 0 and
                math.isfinite(self.degeneracy) and self.degeneracy > 0 and
                v.ndim == 1 and len(v) >= 3 and np.all(np.isfinite(v)) and
                abs(v[0]-self.mass) < 1e-12 and np.all(np.diff(v)>0)):
            raise ValueError('physical ordered energy cell edges required')
        object.__setattr__(self,'edges',v)

    @property
    def size(self):
        return 2*(len(self.edges)-1)


def cell_hats(mode,E):
    z=mode.edges
    E=float(E)
    if not math.isfinite(E) or E<z[0]-1e-12 or E>z[-1]+1e-12:
        raise ValueError('on-shell event outside represented energy cells')
    E=max(z[0],min(z[-1],E))
    j=int(np.searchsorted(z,E,side='right'))-1
    j=max(0,min(j,len(z)-2))
    b=(E-z[j])/(z[j+1]-z[j])
    return 2*j,1-b,b


def mode_thermal_eta(mode,T=1.,mu=0.):
    if T<=0: raise ValueError('positive temperature')
    if not mode.fermion and mu>=mode.mass:
        raise ValueError('Bose chemical potential must remain below mass')
    return ((np.column_stack((mode.edges[:-1],mode.edges[1:])).ravel()-mu)/T)


def mode_moments(mode,eta,order=24):
    """Physical p²dp quantum cell moments, 2x2 block susceptibility."""
    eta=np.asarray(eta,dtype=float)
    if eta.shape!=(mode.size,) or not np.all(np.isfinite(eta)) or order<6:
        raise ValueError('correct finite cell inverse fugacity')
    if not mode.fermion and np.any(eta<=0):
        raise ValueError('strict Bose inverse fugacity for moment inversion')
    xx,ww=leggauss(order)
    M=np.zeros(mode.size)
    J=np.zeros((mode.size,mode.size))
    ent=0.
    fac=mode.degeneracy/(2*math.pi**2)
    for j,(Elo,Ehi) in enumerate(zip(mode.edges[:-1],mode.edges[1:])):
        plo=math.sqrt(max(0,Elo*Elo-mode.mass*mode.mass))
        phi=math.sqrt(max(0,Ehi*Ehi-mode.mass*mode.mass))
        p=(xx+1)*(phi-plo)/2+plo
        E=np.hypot(p,mode.mass)
        b=(E-Elo)/(Ehi-Elo)
        a=1-b
        lam=a*eta[2*j]+b*eta[2*j+1]
        if mode.fermion:
            occ=expit(-lam)
            variance=occ*(1-occ)
            s=np.logaddexp(0,-lam)+lam*occ
        else:
            if np.any(lam<=0): raise ValueError('outside positive Bose domain')
            occ=1/np.expm1(lam)
            variance=occ*(1+occ)
            s=-np.log(-np.expm1(-lam))+lam*occ
        w=fac*ww*p*p*(phi-plo)/2
        M[2*j]+=np.dot(w,a*occ)
        M[2*j+1]+=np.dot(w,b*occ)
        ent+=np.dot(w,s)
        J[2*j,2*j]-=np.dot(w,a*a*variance)
        J[2*j+1,2*j+1]-=np.dot(w,b*b*variance)
        cross=np.dot(w,a*b*variance)
        J[2*j,2*j+1]-=cross
        J[2*j+1,2*j]-=cross
    endpoints=np.column_stack((mode.edges[:-1],mode.edges[1:])).ravel()
    return dict(M=M,J=J,N=float(M.sum()),U=float(np.dot(M,endpoints)),
                S=float(ent),node_energy=endpoints)


def kinematics(p,c,m1,m2,mphi):
    lam=(m1*m1-(m2+mphi)**2)*(m1*m1-(m2-mphi)**2)
    if not(lam>0 and p>=0 and -1<=c<=1):raise ValueError('closed decay')
    e1=math.hypot(p,m1)
    q=math.sqrt(lam)/(2*m1)
    e2=(e1*(m1*m1+m2*m2-mphi*mphi)/(2*m1*m1)+p*q*c/m1)
    eph=e1-e2
    return e1,e2,eph


def scalar_yukawa_width(m1,m2,mphi,y):
    lam=(m1*m1-(m2+mphi)**2)*(m1*m1-(m2-mphi)**2)
    return (y*y*((m1+m2)**2-mphi*mphi)*math.sqrt(lam)/(16*math.pi*m1**3)
            if lam>0 and y>=0 else 0.)


def parent_events(*,m1=.1,m2=.02,mphi=86e-9,y=.001,pmax=2.,nr=12,nangle=12,g1=4):
    if not(pmax>0 and nr>=2 and nangle>=2):raise ValueError
    px,pw=leggauss(nr)
    cs,cw=leggauss(nangle)
    width=scalar_yukawa_width(m1,m2,mphi,y)
    out=[]
    for x,w in zip(px,pw):
        p=pmax*(x+1)/2
        dp=w*pmax/2
        E=math.hypot(p,m1)
        rate=g1/(2*math.pi**2)*dp*p*p*width*m1/E
        out.extend((float(p),float(c),float(rate*w_c/2)) for c,w_c in zip(cs,cw))
    return out,width


def source(modes,eta,events,derivative=False):
    """Same discontinuous hats in F/R, N/E and entropy deposition."""
    if len(modes)!=3 or len(eta)!=3 or not modes[0].fermion or not modes[1].fermion or modes[2].fermion:
        raise ValueError('Fermi Fermi Bose ordering')
    sizes=[m.size for m in modes]
    offsets=np.cumsum([0]+sizes)
    C=np.zeros(sum(sizes))
    jac=np.zeros((len(C),len(C))) if derivative else None
    H=0.
    for p,c,w in events:
        if not(math.isfinite(w) and w>=0):raise ValueError('nonnegative weight')
        en=kinematics(p,c,*[m.mass for m in modes])
        hh=[]; vals=[]
        for m,z,E in zip(modes,eta,en):
            ix,a,b=cell_hats(m,E)
            h=np.zeros(m.size)
            h[ix]=a; h[ix+1]=b
            hh.append(h)
            vals.append(float(np.dot(h,z)))
        f1=expit(-vals[0]); f2=expit(-vals[1])
        if vals[2]<=0: raise ValueError('Bose inverse fugacity nonpositive')
        b=1/math.expm1(vals[2])
        forward=f1*(1-f2)*(1+b)
        reverse=(1-f1)*f2*b
        flux=forward-reverse
        affinity=-vals[0]+vals[1]+vals[2]
        H+=w*flux*affinity
        deriv=(-f1*(1-f1)*(1+b-f2),f2*(1-f2)*(f1+b),-b*(1+b)*(f1-f2))
        for i,sgn in enumerate((-1,1,1)):
            sli=slice(offsets[i],offsets[i+1])
            C[sli]+=w*flux*sgn*hh[i]
            if derivative:
                for j in range(3):
                    slj=slice(offsets[j],offsets[j+1])
                    jac[sli,slj]+=sgn*w*deriv[j]*np.outer(hh[i],hh[j])
    return C,jac,float(H)


def np_block_diag(mats):
    n=sum(len(x) for x in mats)
    out=np.zeros((n,n))
    k=0
    for x in mats:
        m=len(x);out[k:k+m,k:k+m]=x;k+=m
    return out


def ledger(modes,eta,order=24):
    res=[mode_moments(m,z,order) for m,z in zip(modes,eta)]
    return {'M':np.concatenate([x['M'] for x in res]),
            'J':np_block_diag([x['J'] for x in res]),
            'N12':res[0]['N']+res[1]['N'],
            'N1phi':res[0]['N']+res[2]['N'],
            'U':sum(x['U'] for x in res),'S':sum(x['S'] for x in res),
            'Nphi':res[2]['N']}


def implicit_step(modes,eta_old,events,dt,*,order=24,tol=5e-12,max_newton=35):
    if not(dt>0 and math.isfinite(dt)):raise ValueError
    old=ledger(modes,eta_old,order)
    sizes=[m.size for m in modes]
    offsets=np.cumsum([0]+sizes)
    phi_slice=slice(offsets[2],offsets[3])
    state=np.concatenate(eta_old).astype(float)
    if np.any(state[phi_slice]<=0):raise ValueError('Bose state outside physical cone')
    base=old['M']
    scale=max(float(max(abs(base))),1e-12)
    def slices(v):return [v[offsets[i]:offsets[i+1]] for i in range(3)]
    for it in range(max_newton):
        split=slices(state)
        mm=ledger(modes,split,order)
        C,JC,H=source(modes,split,events,True)
        r=mm['M']-base-dt*C
        norm=float(max(abs(r)))
        if norm<tol*scale:
            delta={k:mm[k]-old[k] for k in ('N12','N1phi','U','S','Nphi')}
            lower=float(np.dot(state,mm['M']-base))
            if H < -2e-18 or delta['S'] < -2e-14 or delta['S'] < lower-2e-14:
                raise ArithmeticError('measured entropy or discrete concavity violation')
            return {'eta':split,'iterations':it,'residual':norm,'delta':delta,
                    'instant_H':H,'entropy_lower_bound':lower}
        jac=mm['J']-dt*JC
        direction=solve(jac,-r,assume_a='gen',check_finite=True)
        alpha=1.
        accepted=False
        for _ in range(48):
            cand=state+alpha*direction
            if np.all(np.isfinite(cand)) and np.all(cand[phi_slice]>0):
                try:
                    cs=slices(cand)
                    trial=ledger(modes,cs,order)
                    trialC,_,_=source(modes,cs,events)
                    newnorm=max(abs(trial['M']-base-dt*trialC))
                    if newnorm<norm*(1-1e-4*alpha):
                        state=cand;accepted=True;break
                except (ValueError,ArithmeticError,FloatingPointError):
                    pass
            alpha*=.5
        if not accepted:raise ArithmeticError('Newton cannot descend in positive Bose cone')
    raise ArithmeticError('Newton maximum iteration')


def fixture(*,n=16,mphi=86e-9,pmax=2.,quadratic=True):
    upper=math.hypot(pmax,.1)+.2
    x=np.linspace(0,1,n+1)
    modes=[]
    for mass,fermion in ((.1,True),(.02,True),(mphi,False)):
        edges=mass+(upper-mass)*(x*x if (quadratic and not fermion) else x)
        modes.append(CellMode(mass,edges,fermion=fermion,degeneracy=4 if fermion else 1))
    eta=[mode_thermal_eta(m,.9 if i==0 else 1.) for i,m in enumerate(modes)]
    events,_=parent_events(mphi=mphi,pmax=pmax)
    return modes,eta,events


def tangent_audit(n=16,mphi=86e-9):
    modes,eta,events=fixture(n=n,mphi=mphi)
    C,_,H=source(modes,eta,events)
    starts=np.cumsum([0]+[m.size for m in modes])
    Cp=C[starts[2]:]
    J=mode_moments(modes[2],eta[2],24)['J']
    rate=solve(J,Cp,assume_a='gen')
    return {'n':n,'mphi':mphi,'Cphi_first_cell':Cp[:2].tolist(),
            'eta_rate_first_cell':rate[:2].tolist(),
            'phi_net_rate':float(Cp.sum()),'H':H,
            'tangent_min':float(min(rate)),
            'all_phi_sources_nonpositive':bool(np.all(Cp<=1e-27)),
            'no_real_BEC_claim':True}


if __name__=='__main__':
    import json
    x=tangent_audit()
    print(json.dumps(x,indent=2))
