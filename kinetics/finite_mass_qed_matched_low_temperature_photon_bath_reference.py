"""Low-temperature rest-bath tests with real and virtual photon terms."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
import mpmath as mp

from early_pole_hard_yukawa_photoproduction_reference import G,I,adj,dot,slash
from finite_mass_qed_signed_real_soft_recoil_reference import angular_scaled,setup

L=(I-1j*G[0]@G[1]@G[2]@G[3])/2


def phase(E,m,s):
    return np.sqrt((E*E-(m+s)**2)*(E*E-(m-s)**2))/(2*E)


def trace_born(P,B,scale_m,m,current):
    return np.trace((slash(B)+m*I)@current@(slash(P)+scale_m*I)@adj(current),axis1=-2,axis2=-1).real


def real(m,s,mode,n):
    z,w,p,E,B=setup(1.,m,s,n)
    current=I if mode=='scalar' else L
    P=np.array([1.,0.,0.,0.])
    norm=p*trace_born(P,np.array([E,0.,0.,p]),1.,m,current)
    def f(r):
        S=1-2*r
        sr=np.sqrt(S)
        d=phase(sr,m,s)
        e=np.hypot(d,m)
        gamma=(1-r)/sr
        beta=-r/(1-r)
        b=np.stack([gamma*(e+beta*d*z),d*np.sqrt(1-z*z),0*z,gamma*(d*z+beta*e)],axis=1)
        # k is signed for emission/absorption, with its common r
        # removed from the propagator denominators before squaring.
        kn=np.array([-1.,0.,0.,-1.])
        K=r*slash(kn)
        eik_b=-(m*m/(e+d)+d*(1-z))/sr
        result=0.
        for pol in np.eye(4)[1:3]:
            ep=slash(pol)
            a=current@(2*dot(P,pol)*I+K@ep)/(2*dot(P,kn))
            v=(2*dot(b,pol)[:,None,None]*I-ep@K)@current/(-2*eik_b[:,None,None])
            x=a+v
            result+=np.trace((slash(b)+m*I)@x@(slash(P)+I)@adj(x),axis1=-2,axis2=-1).real
        return (w@result)*d/sr/norm
    if mode=='scalar':
        exact=lambda r:(w@angular_scaled(r,z,1.,m,s,np.sqrt))/norm
        err=max(abs(exact(r)/f(r)-1) for r in (-.01,0,.008))
        if err>3e-8:raise ArithmeticError(err)
    return f,B


def shifted_rate(m,s,mode):
    # h = (g²/4π²) J1 T² in M=1 units. Both pole mass squares gain 4h.
    # Eq. (24),(30),(33) of Czarnecki et al.: use the spin vector p-K
    # and the unshifted mass in the on-shell spinor numerator.
    with mp.workdps(55):
        m,s=map(mp.mpf,(m,s))
        def f(h,part):
            M1,m1=mp.sqrt(1+4*h),mp.sqrt(m*m+4*h)
            E=(M1*M1+m1*m1-s*s)/(2*M1)
            p=mp.sqrt((E-m1)*(E+m1))
            K0p=2*h/M1
            K0=h*mp.log((E+p)/(E-p))/p
            if part=='phase':
                # Isolate the phase-space p/M_T²: freeze the Born trace.
                return p/M1**2
            if part=='trace':
                return (M1-K0p)*(E-K0)+(m if mode=='scalar' else 0)
            return p/M1**2*((M1-K0p)*(E-K0)+(m if mode=='scalar' else 0))
        values={}
        for part in ('phase','trace','combined'):
            values[part]=float(mp.diff(lambda h:f(h,part),0)/f(0,part))
        return values



def virtual_eikonal(m,s,mode):
    M=1.
    p=phase(M,m,s)
    E=np.hypot(p,m)
    P=np.array([1.,0.,0.,0.])
    B=np.array([E,0.,0.,p])
    current=I if mode=='scalar' else L
    U,V=slash(B)+m*I,slash(P)+I
    reference=trace_born(P,B,1.,m,current)
    d=dot(P,B)
    worst=0.
    for w in (1e-5,.002,.03):
        for u in (-.81,-.17,.31,.94,.9992):
            K=w*np.array([1.,np.sqrt(1-u*u),0.,u])
            a,b=dot(P,K),dot(B,K)
            vertex=0.
            for sign in (-1,1):
                J=slash(sign*K)
                Gamma=sum(g*Z@(U-J)@current@(V-J)@Z for g,Z in zip((1.,-1.,-1.,-1.),G,strict=True))
                vertex+=np.trace(U@Gamma@V@adj(current)).real/(reference*4*a*b)
            real=2*(2*d/(a*b)-M*M/a**2-m*m/b**2)
            pole=2*(M*M/a**2+m*m/b**2)
            physical_vertex=-2*vertex
            err=abs(real+pole+physical_vertex)/max(1.,abs(real),abs(pole),abs(physical_vertex))
            worst=max(worst,float(err))
    return worst


def c2_scalar_complex(m,s,n):
    z,w,p,E,_=setup(1.,m,s,n)
    norm=p*4*(E+m)
    theta=2*np.pi*np.arange(40)/40
    out=[]
    for r in (.003,.006):
        xs=r*np.exp(1j*theta)
        yy=np.array([(w@angular_scaled(x,z,1.,m,s,np.sqrt))/norm for x in xs])
        out.append(2*np.mean((yy-yy[0])*np.exp(-2j*theta)).real/r**2)
    return out


def finite_t(calculation,m,s,mode,ext,tau):
    f,soft=calculation
    c2=ext['real_coefficient']
    kmax=(1-(m+s)**2)/2
    def integrand(t):
        x=t*tau
        if x<4e-5:
            residual=c2*x*x
        else:
            residual=f(-x)+(f(x) if x<kmax else 0.)-2*soft
        return residual/(t*np.expm1(t)) if t>0 else 0.
    cuts=sorted({0.,1.,5.,10.,25.,60.,4e-5/tau,*([kmax/tau] if kmax/tau<60 else [])})
    thermal=0.
    for a,b in zip(cuts[:-1],cuts[1:],strict=True):
        thermal+=quad(integrand,a,b,epsabs=3e-13,epsrel=1e-7,limit=100)[0]
    J1=np.pi**2/6
    return dict(temperature_over_parent_mass=tau,raw_real_subtracted=thermal,
                normalized_full_thermal_correction=thermal/(J1*tau*tau)+ext['external']['combined'])


def calc(m,s,mode,do_finite=False):
    f,B=real(m,s,mode,84)
    b0=f(0.)
    if abs(b0-B)>2e-8:
        raise ArithmeticError('Real matrix element missed the soft limit.')
    def ladder(fun,steps):
        prev=None
        results=[]
        soft=fun(0.)
        for h in steps:
            fd=(fun(h)+fun(-h)-2*soft)/(h*h)
            rich=(4*fd-prev)/3 if prev is not None else None
            results.append(dict(step=h,even_second_derivative=float(fd),richardson=float(rich) if rich is not None else None))
            prev=fd
        return results
    first=ladder(f,(.002,.001,.0005))
    other,_=real(m,s,mode,116)
    second=ladder(other,(.0016,.0008,.0004))
    a=first[-1]['richardson']
    b=second[-1]['richardson']
    complex_replay=c2_scalar_complex(m,s,100) if mode=='scalar' else None
    if abs(a-b)>1.5e-6 or (complex_replay is not None and max(abs(v-b) for v in complex_replay)>1.5e-6):
        raise ArithmeticError('Low-T recoil coefficients did not independently reproduce.')
    ext=shifted_rate(m,s,mode)
    if abs(ext['phase']+ext['trace']-ext['combined'])>1e-12:
        raise ArithmeticError('External phase and spin coefficients did not reconstruct the combined derivative.')
    result=a+ext['combined']
    check=virtual_eikonal(m,s,mode)
    if check>2e-8 or abs(result+2)>1.5e-6:
        raise ArithmeticError(f'Matched thermal equal-charge coefficient failed: {result}, {check}.')
    out=dict(current=mode,m_over_M=m,scalar_over_M=s,soft_angular=B,
             real_coefficient=a,independent_derivative=b,scalar_complex_contour=complex_replay,
             first_ladder=first,independent_ladder=second,external=ext,total_T2_in_J1_units=result,
             photon_real_vertex_pole_soft_defect=check)
    if do_finite:
        ts=(.001,.002,.004)
        out['finite_temperature_ladder']=[finite_t((f,b0),m,s,mode,out,tau) for tau in ts]
        errors=[abs(v['normalized_full_thermal_correction']-result) for v in out['finite_temperature_ladder']]
        if max(errors)>.02 or not errors[0]<errors[1]<errors[2]:
            raise ArithmeticError(f'Finite-temperature residual did not converge down the bath ladder: {errors}.')
    return out


def run():
    cases=[]
    for current in ('left','scalar'):
        for m,s in ((.1,0.),(.03,0.),(.01,0.),(.003,0.),(.2,.1),(.5,.12),(.06,.4),(.2,.65)):
            cases.append(calc(m,s,current))
    archive=[]
    for name in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{name}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        archive.append(dict(grid=name,**calc(d['mass_ratio'],d['checkpoints'][-1]['scalar_pole_over_T']/d['x_end'],
                                            'scalar',True)))
    return dict(scope='Isotropic low-T photon-only rest-bath control, equal charges and no thermal final-state factors or charged-particle occupations. Not a high-T archive application.',
                normalization='multiply the T^2 coefficient by g^2/(4 pi^2) (pi^2/6) (T/M)^2 relative to Born',
                archive_masses_only=archive,independent_mass_and_current_controls=cases)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
