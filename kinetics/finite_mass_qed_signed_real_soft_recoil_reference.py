"""Signed exact real recoil plus leading eikonal subtraction for a rest-frame control."""
import gzip
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import amplitudes


def setup(M,m,s,n):
    p0=np.sqrt((M*M-(m+s)**2)*(M*M-(m-s)**2))/(2*M)
    Eb=np.hypot(p0,m)
    beta=p0/Eb
    l,r=np.log((m/Eb)**2/(1+beta)),np.log1p(beta)
    z,w=roots_legendre(n)
    u=(l+r)/2+(r-l)*z/2
    x=(1-np.exp(u))/beta
    aw=w*(r-l)/4*np.exp(u)/beta
    return np.asarray(x),np.asarray(aw),p0,Eb,2*(np.arcsinh(p0/m)/beta-1)


def angular_scaled(r,x,M,m,s,sqrt):
    k=M*r
    S=M*(M-2*k)
    root=sqrt(S)
    p=sqrt((S-(m+s)**2)*(S-(m-s)**2))/(2*root)
    E=sqrt(p*p+m*m)
    beta=-k/(M-k)
    gamma=(M-k)/root
    Y=gamma*(E+beta*p*x)
    Z=gamma*(p*x+beta*E)
    A,C=-M,-(Y-Z)
    d=M*Y
    # p_gamma^2 |M_real|^2 after canceling the singular denominators
    # algebraically; all factors here have a smooth signed r -> 0 limit.
    V=(1/A-1/C)/2
    H=M*M/(A*A)+m*m/(C*C)-2*d/(A*C)
    X=C*M*M/A+A*m*m/C-2*d
    U=-4*(d+M*m)*H-8*k*V*X+16*k*k*V*V*A*C
    return (p/root)*U


def trace_check(M,m,s):
    out=[]
    z=np.asarray([-.85,-.2,.3,.74,.96])
    for r in (-.01,-1e-4,1e-4,.01):
        k=M*r
        ss=M*(M-2*k)
        w=np.sqrt(ss)
        p=np.sqrt((ss-(m+s)**2)*(ss-(m-s)**2))/(2*w)
        E=np.hypot(p,m)
        ga=(M-k)/w
        be=-k/(M-k)
        b=np.array([ga*(E+be*p*z),p*np.sqrt(1-z*z),0*z,ga*(p*z+be*E)]).T
        incoming=np.array([M,0.,0.,0.])
        signed=np.array([-k,0.,0.,-k])
        value=sum(amplitudes(incoming,signed,b,M,m,u)[0] for u in np.eye(4)[1:3])
        got=angular_scaled(r,z,M,m,s,np.sqrt)
        out.append(float(max(abs(got/(value*k*k*p/w)-1))))
    if max(out)>2e-10:
        raise ValueError('Signed external photon spin-trace check failed.')
    return max(out)


def sample(M,m,s,n):
    x,w,p0,E0,B=setup(M,m,s,n)
    norm=(p0/M)*4*M*(E0+m)
    def direct(r):
        return float((w@angular_scaled(r,x,M,m,s,np.sqrt))/norm)
    with mp.workdps(65):
        ma,mb,ms=map(mp.mpf,(M,m,s))
        ab=list(zip(map(mp.mpf,x.tolist()),map(mp.mpf,w.tolist()),strict=True))
        pp=mp.sqrt((ma*ma-(mb+ms)**2)*(ma*ma-(mb-ms)**2))/(2*ma)
        nn=(pp/ma)*4*ma*(mp.sqrt(mb*mb+pp*pp)+mb)
        def exact(r):
            return mp.fsum(w*angular_scaled(r,x,ma,mb,ms,mp.sqrt) for x,w in ab)/nn
        c=[float(mp.diff(exact,mp.mpf(0),j)/mp.factorial(j)) for j in range(5)]
        verify=[]
        for r in (mp.mpf('1e-4'),mp.mpf('1e-7')):
            exact_odd=exact(r)-exact(0)
            exact_even=exact(r)+exact(-r)-2*exact(0)
            odd=sum(mp.mpf(c[j])*r**j for j in range(1,5))
            even=2*sum(mp.mpf(c[j])*r**j for j in (2,4))
            verify.append([float(abs(odd/exact_odd-1)),float(abs(even/exact_even-1))])
    if abs(c[0]-B)>1e-10 or max(max(v) for v in verify)>5e-10:
        raise ValueError('The signed real recoil did not recover the eikonal or expansion test.')

    def rates(r):
        bn=1/np.expm1(M*r/.9)
        if r<=1e-4:
            one=sum(c[j]*r**j for j in range(1,5))
            two=2*(c[2]*r*r+c[4]*r**4)
        else:
            a,b=direct(r),direct(-r)
            one=a-c[0]
            two=a+b-2*c[0]
        return one/r,bn*two/r

    def fw(r,k):return rates(r)[k]
    comps=[]
    for cut in (0.,1e-4,1e-6,1e-10):
        v=[quad(lambda r:fw(r,k),cut,.01,points=[1e-4] if cut<1e-4 else None,
                epsabs=3e-10,epsrel=3e-10)[0] for k in (0,1)]
        comps.append(dict(normalized_ir_cut=cut,vacuum_residual=v[0],thermal_even_residual=v[1],
                          sum_after_eikonal_subtraction=float(sum(v))))
    return dict(normalized_soft_coefficients=c,analytic_zero_soft=B,
                max_series_test_relative_defect=max(max(v) for v in verify),
                angle_nodes=n,cutoff_ladder=comps)


def run():
    rows=[]
    for grid in ('fine','third'):
        path=Path(__file__).parent/f'results/early-pole-spectrum-{grid}-m06-2026-11-26.json.gz'
        with gzip.open(path) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        a=sample(M,m,s,28)
        b=sample(M,m,s,52)
        if abs(a['cutoff_ladder'][0]['sum_after_eikonal_subtraction']-
               b['cutoff_ladder'][0]['sum_after_eikonal_subtraction'])>5e-8:
            raise ValueError('Collinear angular quadrature did not resolve the paired signed residual.')
        rows.append(dict(grid=grid,independent_signed_trace_defect=trace_check(M,m,s),coarse=a,fine=b))
    return dict(scope='Exact signed vacuum matrix recoil with a parent-rest photon bath and leading-soft virtual counterterm only.',
                upper_photon_energy_over_parent_mass=.01,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
