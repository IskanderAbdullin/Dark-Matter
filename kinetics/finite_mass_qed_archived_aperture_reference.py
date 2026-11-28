"""Frozen, unweighted parent-momentum aperture with a real O(g²) dispersion."""
import functools
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import intervals, Q


class Frozen:
    def __init__(self,d,name,m,cut):
        self.edges=np.asarray(d[f'endpoint_{name}_momentum_edges'])
        self.f=np.asarray(d[f'endpoint_{name}_occupations'])
        self.m=m
        v=self.edges**3/3
        self.center=(v[1:]+v[:-1])/2
        sec=np.diff(self.f)/np.diff(self.center)
        slope=np.zeros_like(self.f)
        slope[1:-1]=np.where(sec[:-1]*sec[1:]>0,
                            np.sign(sec[:-1])*np.minimum(abs(sec[:-1]),abs(sec[1:])),0.)
        limit=2*np.minimum(self.f,1-self.f)/np.diff(v)
        self.slope=np.clip(slope,-limit,limit)
        i=np.searchsorted(self.edges,cut,side='right')-1
        edges=self.edges[i:].copy()
        edges[0]=cut
        self.lo,self.hi,self.fhard=edges[:-1],edges[1:],self.f[i:]
        self.photon=np.pi**2/6-quad(lambda q:q/np.expm1(q),0.,cut,epsabs=1e-13)[0]
        t,w=roots_legendre(24)
        self.r=(self.lo[:,None]+self.hi[:,None])/2+np.diff(edges)[:,None]*t/2
        v=np.hypot(self.r,m)
        self.w=np.diff(edges)[:,None]*w/2*self.fhard[:,None]
        self.moment=float(np.sum(self.w*self.r**2/v))
        self.amp=self.w*self.r/v

    def f_at(self,p,limited=False):
        j=np.searchsorted(self.edges,p,side='right')-1
        if j<0 or j>=len(self.f):
            raise ValueError('Aperture escaped input support.')
        return self.f[j]+self.slope[j]*(p**3/3-self.center[j]) if limited else self.f[j]

    @functools.lru_cache(maxsize=8192)
    def log(self,p):
        # Integrate the log exactly to quadrature precision around p. Far cells
        # are smooth; close cells use a direct adaptive or singular-weight rule.
        width=self.hi-self.lo
        near=(self.lo-2*width <= p)&(p <= self.hi+2*width)
        gap=abs(self.r[~near]-p)
        total=np.sum(self.amp[~near]*np.log1p(2*np.minimum(self.r[~near],p)/gap))
        for l,h,f in zip(self.lo[near],self.hi[near],self.fhard[near],strict=True):
            def a(z): return z/np.hypot(z,self.m)
            if l<=p<=h:
                val=quad(lambda z:a(z)*np.log(z+p),l,h,epsabs=2e-13)[0]
                for sign,L in ((-1,p-l),(1,h-p)):
                    if L>0:
                        val-=quad(lambda s:a(p+sign*s),0,L,weight="alg-loga",
                                  wvar=(0.,0.),epsabs=2e-13)[0]
            else:
                val=quad(lambda z:a(z)*np.log1p(2*min(z,p)/abs(z-p)),l,h,
                         epsabs=2e-13)[0]
            total+=f*val
        return -self.m**2*total/(2*p)

    def factor(self,p,logs):
        a=self.photon+self.moment+(self.log(p) if logs else 0.)
        return a/(2*np.pi**2*np.hypot(p,self.m))

    def energy(self,p,g2,logs):
        return np.hypot(p,self.m)+g2*self.factor(p,logs)


def difference(P,p,M,b,Es):
    return ((M-b)*(M+b)+(P-p)*(P+p))/(np.hypot(P,M)+np.hypot(p,b))-Es


def solve_window(parent,daughter,ms,Q,g2,logs):
    M,b=parent.m,daughter.m
    Es,E0,E1,*_=intervals(M,ms,b,Q)
    P0,P1=np.sqrt((np.asarray((E0,E1))-M)*(np.asarray((E0,E1))+M))
    # Separate collinear residue from subtraction of two large pole energies.
    def unfold(P,k):
        return difference(P,k,M,b,Es)+g2*(parent.factor(P,logs)-daughter.factor(k,logs))
    def face_res(P,s):
        # Compute P-k exactly on the collinear faces rather than losing digits.
        k=P+s*Q
        r=((M-b)*(M+b)-2*s*P*Q-Q*Q)/(np.hypot(P,M)+np.hypot(k,b))-Es
        return r+g2*(parent.factor(P,logs)-daughter.factor(k,logs))
    def root(fn,seed,step=.018):
        return brentq(fn,seed-step,seed+step,xtol=9e-15,rtol=1e-15)
    Plo,Phi=(root(lambda P:face_res(P,s),Pc) for s,Pc in ((1,P0),(-1,P1)))
    if g2==0 and max(abs(Plo-P0),abs(Phi-P1))>1e-12:
        raise ValueError('Vacuum roots disagree with the closed two-mass formula.')
    if not Plo<Phi or max(abs(Plo-P0),abs(Phi-P1))>.012:
        raise ValueError('Wrong aperture root branch.')
    slopecheck=np.linspace(Plo-1e-4,Phi+1e-4,15)
    if any(np.any(np.diff([face_res(z,s) for z in slopecheck])>=0) for s in (-1,1)):
        raise ValueError('Collinear hard residues are not monotonically decreasing.')
    if Plo-Q <= max(.25,10*np.sqrt(g2)):
        raise ValueError('Daughter support leaked outside the tested hard domain.')
    rhs=[np.asarray([F.energy(z,g2,logs) for z in np.linspace(Plo-Q,Phi+Q,29)])
         for F in (parent,daughter)]
    if any(np.any(np.diff(u)<=0) for u in rhs):
        raise ValueError('Energy inversion is not monotone in the hard interval.')
    # Daughter shell edges push through the energy-matching map, not through
    # a constant vacuum-mass or common-momentum approximation.
    breaks=[Plo,Phi]
    breaks+=parent.edges[(parent.edges>Plo)&(parent.edges<Phi)].tolist()
    for k in daughter.edges[(daughter.edges>=Plo-Q)&(daughter.edges<=Phi+Q)]:
        if unfold(Plo,k)<0<unfold(Phi,k):
            breaks.append(brentq(lambda z:unfold(z,k),Plo,Phi,xtol=9e-15,rtol=1e-15))
    br=np.unique(breaks)
    seed=stim=0.
    smooth=np.zeros((2,2))
    minstim=np.inf
    defect=0.
    for left,right in zip(br[:-1],br[1:],strict=True):
        vals=[]
        for p in (left+0.13*(right-left),(left+right)/2,right-0.13*(right-left)):
            k=brentq(lambda k:unfold(p,k),p-Q,p+Q,xtol=9e-15,rtol=1e-15)
            defect=max(defect,abs(unfold(p,k)))
            fp,fb=parent.f_at(p),daughter.f_at(k)
            vals.append((fp,fb))
        if len(set(vals))!=1:
            raise ValueError('An occupation cell crossing was missed.')
        fp,fb=vals[0]
        seed+=(right-left)*fp*(1-fb)
        stim+=(right-left)*(fp-fb)
        minstim=min(minstim,float(fp-fb))
        # The same exact occupation-cell cuts apply to finite-volume
        # minmod polynomials in u=p³/3 which retain shell number.
        for i,nq in enumerate((4,8)):
            t,w=roots_legendre(nq)
            for p,weight in zip((left+right)/2+(right-left)*t/2,w,strict=True):
                k=brentq(lambda k:unfold(p,k),p-Q,p+Q,xtol=9e-15,rtol=1e-15)
                fp,fb=parent.f_at(p,True),daughter.f_at(k,True)
                if not (0<=fp<=1 and 0<=fb<=1):
                    raise ValueError('Limited occupancy inadmissible.')
                smooth[i]+=weight*(right-left)/2*np.asarray([fp*(1-fb),fp-fb])
    if max(abs(smooth[0]-smooth[1]))>2e-13:
        raise ValueError('Limited reconstruction mesh integration mismatch.')
    if defect>1e-14:
        raise ValueError('Hard daughter energy matching failed.')
    return dict(endpoints=[Plo,Phi],center=(Plo+Phi)/2,width=Phi-Plo,splits=len(br)-2,
                unweighted_seed=seed,unweighted_stimulus=stim,
                minimum_pointwise_stimulus=minstim,root_max_energy_defect=defect,
                limited_seed=float(smooth[1,0]),limited_stimulus=float(smooth[1,1]),
                limited_quadrature_defect=float(max(abs(smooth[0]-smooth[1]))))


def run():
    out=[]
    inputs={s:json.load(gzip.open(Path(__file__).parent/
                                f'results/early-pole-spectrum-{s}-m06-2026-11-26.json.gz'))
            for s in ('fine','third')}
    poles={s:d['checkpoints'][-1]['scalar_pole_over_T'] for s,d in inputs.items()}
    for name,d in inputs.items():
        M=d['x_end'];b=M*d['mass_ratio']
        for cut in (.1,.25,.5):
            parent,daughter=Frozen(d,'parent',M,cut),Frozen(d,'fermion',b,cut)
            for pole,ms in poles.items():
                bare=solve_window(parent,daughter,ms,Q,0.,False)
                slopes=[]
                for s,P in zip((1,-1),bare['endpoints'],strict=True):
                    k=P+s*Q
                    ep,eb=np.hypot(P,M),np.hypot(k,b)
                    rb=(-2*s*Q-np.hypot(ms,Q)*(P/ep+k/eb))/(ep+eb)
                    slopes.append(-(parent.log(P)/(2*np.pi**2*ep)-daughter.log(k)/(2*np.pi**2*eb))/rb)
                variants=[]
                for g in (.001,.005,.01):
                    fac=solve_window(parent,daughter,ms,Q,g*g,False)
                    full=solve_window(parent,daughter,ms,Q,g*g,True)
                    relative_shift=np.asarray(full['endpoints'])-fac['endpoints']
                    slope_error=np.max(abs(relative_shift/(g*g*np.asarray(slopes))-1))
                    if g==.001 and slope_error>.001:
                        raise ValueError('Nonlinear validation of bare log slope failed.')
                    # Independent comparison for the constant-mass-squared
                    # model. The real O(g²) energy expansion differs from its
                    # exact resummation by O(g⁴), so use g-scaled defects.
                    md,bm=[np.sqrt(F.m*F.m+g*g*(F.photon+F.moment)/np.pi**2)
                           for F in (parent,daughter)]
                    _,l,h,*_=intervals(md,ms,bm,Q)
                    exact=np.sqrt((np.asarray([l,h])-md)*(np.asarray([l,h])+md))
                    mass_check=max(abs(exact-np.asarray(fac['endpoints'])))
                    if mass_check>=2e3*g**4:
                        raise ValueError('Factorized support failed the vacuum-invariant control.')
                    vals={k:(full[k]/fac[k]-1) for k in ('width','unweighted_seed','unweighted_stimulus','limited_seed','limited_stimulus')}
                    variants.append(dict(g=g,factorized=fac,with_log=full,
                                         log_changes_in_bare_width_units=(relative_shift/bare['width']).tolist(),
                                         log_relative_to_factorized=vals,slope_relative_defect=float(slope_error),
                                         factorized_vs_exact_squared_mass_root_defect=float(mass_check)))
                out.append(dict(archive=name,scalar_pole_grid=pole,cut=cut,bare=bare,
                                log_momentum_derivatives_at_g2_zero=slopes,variants=variants))
    return dict(scope='Frozen collinear support and unweighted p integrals, not a physical decay rate.',
                scalar_pivot=Q,rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
