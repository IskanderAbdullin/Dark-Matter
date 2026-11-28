"""Collisionless isotropic transverse photon pole from two frozen charged spectra."""
import gzip
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

Q=4.1e-6
G=.005
PREF=2*G*G/np.pi**2  # two spins and two charge signs PER species, equal anti proxy


@lru_cache(None)
def gl(n):
    return roots_legendre(n)


def transverse(x):
    x=np.asarray(x,dtype=float)
    h=np.zeros_like(x)
    d=np.zeros_like(x)
    low=x<.28
    y=x[low]
    for n in range(15):
        c=1./((2*n+1)*(2*n+3))
        h[low]+=c*y**(2*n)
        if n:d[low]+=2*n*c*y**(2*n-1)
    y=x[~low]
    ah=np.arctanh(y)
    h[~low]=.5*(1/y**2+(1-1/y**2)*ah/y)
    d[~low]=.5*(-3/y**3+(3/y**4-1/y**2)*ah)
    return h,d


def angular_fixtures():
    u,w=gl(240)
    x=np.array([0.,.03,.18,.27,.2801,.4,.8,.95,.990,.997])
    h,d=transverse(x)
    h0=((1-u*u)/(1-x[:,None]*u))@w/4
    d0=((1-u*u)*u/(1-x[:,None]*u)**2)@w/4
    err=float(max(np.max(abs(h0-h)),np.max(abs(d0-d))))
    if err>4e-9:raise ArithmeticError(('Independent polar integral',err))
    # Equilibrium massless spin+charge multiplicity, then finite-mass
    # susceptibility against direct (minus d f / dE) integrals.
    fd=lambda p:expit(-p)*expit(p)
    D=PREF*quad(lambda p:p*p*fd(p),0,70,epsabs=1e-12)[0]
    z=abs(D/(G*G/3)-1.)
    if z>5e-14:raise ArithmeticError('One massless Dirac degree-of-freedom normalization')
    finite=[]
    for M in (.003,.05,.2,1.):
        for T in (.4,1.):
            def eps(x):return np.hypot(x,M)
            a=PREF*quad(lambda p:(p**4/(3*eps(p)**2))*expit(-eps(p)/T)*expit(eps(p)/T)/T,
                       0.,70*T,epsabs=1e-19,epsrel=1e-12)[0]
            b=PREF*quad(lambda p:p*p/eps(p)*(1-p*p/(3*eps(p)**2))*expit(-eps(p)/T),
                       0.,70*T,epsabs=1e-19,epsrel=1e-12)[0]
            finite.append(abs(a/b-1.))
    if max(finite)>2e-10:raise ArithmeticError('Massive plasma-frequency integral identity')
    return dict(angular_value_and_derivative_max_abs_defect=err,massless_Debye_theory_relative_error=z,
                finite_mass_equilibrium_plasma_frequency_identity_relative_defect=max(finite))


class Species:
    def __init__(self,m,edges,f):
        self.m=m
        self.e=np.asarray(edges)
        self.f=np.asarray(f)
        sign=np.diff(np.r_[0.,self.f,0.])
        self.p=self.e
        self.E=np.hypot(self.p,m)
        self.v=self.p/self.E
        self.W=-sign*self.p*self.E
        self.D=PREF*sum(self.W)
        self.G=PREF*sum(self.W*self.v**2)/3

    def field(self,beta):
        h,d=transverse(beta*self.v)
        return PREF*(self.W*self.v**2@h),PREF*(self.W*self.v**3@d)

    def by_parts(self,beta,order):
        x,w=gl(order)
        l,h=self.e[:-1,None],self.e[1:,None]
        p=(l+h)/2+(h-l)*x/2
        E=np.hypot(p,self.m)
        v=p/E
        vp=self.m**2/E**3
        H,D=transverse(beta*v)
        der=(2*p*v+p*p*vp)*H+p*p*v*D*beta*vp
        z=PREF*(self.f*((h-l)[:,0]/2))@(der@w)
        dd=PREF*(self.f*((h-l)[:,0]/2))@((E+p*p/E)@w)
        return z,dd


def case(d,anti_proxy=True):
    M=d['x_end'];m=M*d['mass_ratio']
    s=d['checkpoints'][-1]['scalar_pole_over_T']
    Es=np.hypot(s,Q)
    charges=[Species(a,d[f'endpoint_{key}_momentum_edges'],d[f'endpoint_{key}_occupations'])
             for a,key in ((M,'parent'),(m,'fermion'))]
    fac=1. if anti_proxy else .5
    massD=fac*sum(c.D for c in charges)
    gap2=fac*sum(c.G for c in charges)
    gap=np.sqrt(gap2)
    def pi(beta):
        return fac*np.sum([c.field(beta) for c in charges],axis=0)
    def pole(k):
        if k==0.:return gap,0.,1.
        def z(w):
            return (w-k)*(w+k)-pi(k/w)[0]
        w=brentq(z,max(np.nextafter(k,np.inf),gap),np.sqrt(k*k+massD),xtol=5e-17)
        P,beta_der=pi(k/w)
        dOmega2=-beta_der*k/(2*w**3)
        res=1./(1-dOmega2)
        def dv(w):
            return (w-k)*(w+k)-pi(k/w)[0]
        if abs(dv(w))>5e-13*gap2 or not (0.<res<=1.):
            raise ArithmeticError(('Transverse pole or residue not physical in this fixture',w,res,dv(w)))
        return w,(w-k)*(w+k),res
    if gap<=Es:
        raise ArithmeticError('This archive/coupling fixture does not have the assumed past-timelike start')
    delta=gap*5e-4
    z=pole(delta)[0]
    c2=1+fac*PREF*sum(c.W@c.v**4 for c in charges)/(15*gap2)
    observed=(z-gap)*(z+gap)/delta**2
    if abs(observed/c2-1.)>1e-6:
        raise ArithmeticError('Independent small-momentum plasma pole expansion')
    rows=[]
    max_parts=max_partsD=0.
    for k in (0.,4e-5,.0002,.0004,.0008,.00128,.00256,.005,.01,.02):
        w,inv,R=pole(k)
        iso=w is not None and k==0
        if iso:inv=gap2
        cur=pi(k/w)[0]
        ibp=fac*sum(ch.by_parts(k/w,40)[0] for ch in charges)
        max_parts=max(max_parts,abs(cur/ibp-1.))
        for c in charges:
            max_partsD=max(max_partsD,abs(c.D/c.by_parts(k/w,40)[1]-1.))
        tf=(Es+w)**2-(k+Q)**2
        ta=(Es-w)**2-(k+Q)**2
        tb=(Es-w)**2-(k-Q)**2
        rows.append(dict(k=k,omega=w,omega2_minus_k2=inv,physical_transverse_w2_pole_residue=R,
                         lowest_emission_transfer_invariant=tf,
                         min_max_absorption_transfer_invariant_over_angles=[ta,tb]))
    threshold=[]
    finite_support=[]
    Emax=[float(ch.E[-1]) for ch in charges]
    for c in (-1.,-.5,0.,.5,1.):
        def V(k):
            return np.hypot(k-Q*c,Q*np.sqrt((1-c)*(1+c)))
        def difference(k):
            return pole(k)[0]-Es-V(k)
        turn=brentq(difference,0.,.025,xtol=5e-16)
        if difference(.999*turn)<=0 or difference(1.001*turn)>=0:
            raise ArithmeticError('Absorption past-timelike crossing not bracketed')
        threshold.append([c,turn])
        def parent_min(k):
            e=Es-pole(k)[0]
            u=(e-V(k))*(e+V(k))
            if e>=0. or u>=0.:raise ArithmeticError('Wrong transfer branch')
            D=M*M+u-m*m
            H=np.sqrt(D*D-4*M*M*u)
            # The past null boundary is different from the future null
            # vacuum chart. This expression grows like 1/(-u), not finite.
            return (e*D-V(k)*H)/(2*u)
        edges=[]
        for side in (0,1):
            def bounds(k):
                gain_cap,loss_cap=Emax[0],Emax[1]+Es-pole(k)[0]
                return parent_min(k)-(gain_cap if side==0 else loss_cap)
            ar=brentq(bounds,turn*1.0000001,.025,xtol=6e-16)
            if not ar>turn or abs(bounds(ar))>2e-8:
                raise ArithmeticError('Finite archive support gate')
            edges.append(ar)
        finite_support.append(dict(cosine=c,transfer_spacelike_first_k=turn,
                corresponding_photon_energy=pole(turn)[0],
                first_possible_parent_and_daughter_occupied_k=edges,
                first_possible_parent_and_daughter_occupied_photon_energy=[pole(z)[0] for z in edges],
                required_parent_energy_at_1p001_turn=parent_min(turn*1.001)))
    err=max(max_parts,max_partsD)
    if err>2e-8:raise ArithmeticError(('Exact jump vs radial integration by parts',err))
    return dict(coupling=G,single_sign_if_no_anti_proxy=not anti_proxy,held_scalar_E=Es,
                charged_only_Debye_mass=float(np.sqrt(massD)),plasma_frequency=gap,
                plasma_frequency_over_scalar_energy=gap/Es,
                small_k_squared_slope=c2,small_k_slope_check_relative_defect=abs(observed/c2-1.),
                charge_components_Wpole_squared=[fac*c.G for c in charges],
                derivative_measure_vs_positive_radial_relative_defect=err,
                absorption_spacelike_entry_by_cosine=threshold,finite_archive_support_gate=finite_support,poles=rows)


def run():
    root=Path(__file__).parent/'results'
    rows=[]
    for g in ('fine','third'):
        with gzip.open(root/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        rows.append(dict(grid=g,default=case(d),no_anti_proxy_control=case(d,False)))
    return dict(scope='Collisionless isotropic kinetic polarization of bare massive charged constituents, one Abelian photon, no phenomenological width or projected scattering rate.',fixtures=angular_fixtures(),data=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
