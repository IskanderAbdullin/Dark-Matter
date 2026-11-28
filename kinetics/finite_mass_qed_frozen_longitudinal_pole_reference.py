"""Independent collisionless longitudinal pole and geometry gates from frozen charged occupancies."""
import gzip,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q,PREF,Species,gl


def L(x):
    h=np.empty_like(x);d=np.empty_like(x)
    a=x<.25;y=x[a]
    h[a]=sum(y**(2*j)/(2*j+1) for j in range(1,14))
    d[a]=sum(2*j*y**(2*j-1)/(2*j+1) for j in range(1,14))
    y=x[~a];y0=np.arctanh(y)
    h[~a]=y0/y-1.
    d[~a]=(y/(1-y*y)-y0)/(y*y)
    return h,d


def evaluate(g,data):
    M=data['x_end'];m=M*data['mass_ratio']
    s=data['checkpoints'][-1]['scalar_pole_over_T']
    E=np.hypot(s,Q)
    ch=[Species(v,data[f'endpoint_{lab}_momentum_edges'],data[f'endpoint_{lab}_occupations'])
        for v,lab in ((M,'parent'),(m,'fermion'))]
    wp2=sum(x.G for x in ch)
    wp=np.sqrt(wp2)
    def response(beta):
        return PREF*sum(x.W@L(beta*x.v)[0] for x in ch)
    kL=np.sqrt(response(1.))
    t,w=gl(48)
    def bpart(beta):
        integral=0.
        for c in ch:
            a,b=c.e[:-1,None],c.e[1:,None]
            p=(a+b)/2+(b-a)*t/2
            E=np.hypot(p,c.m);v=p/E
            f,d=L(beta*v)
            dx=(E+p*p/E)*f+p*E*d*beta*c.m**2/E**3
            integral+=PREF*((c.f*(b-a)[:,0]/2)@(dx@w))
        return integral
    cuterr=max(abs(response(x)/bpart(x)-1.) for x in (.04,.22,.6,.92,1.))
    if cuterr>5e-7:raise ArithmeticError(('Longitudinal jump/radial mismatch',cuterr))
    def pole(k):
        if k==0.:return wp
        if k>=kL:raise ValueError('No timelike longitudinal root above light-cone endpoint in this frozen model')
        upper=max(.01,1.2*k)
        for _ in range(10):
            if response(k/upper)<k*k:break
            upper*=2
        return brentq(lambda w:response(k/w)-k*k,np.nextafter(k,np.inf),upper,xtol=5e-17)
    d=wp*7e-4
    z=pole(d)
    slope=PREF*sum(c.W@c.v**4 for c in ch)/(5*wp2)
    accepted=(z-wp)*(z+wp)/d**2
    if abs(accepted/slope-1.)>2e-6:raise ArithmeticError('Small-k slope mismatch')
    kk=[0.,.0002,.0004,.0008,.00128,.00256,.005]
    kk=sorted(set([x for x in kk if x<kL]+[.999*kL]))
    rows=[]
    for q in kk:
        w=pole(q)
        rows.append(dict(k=q,omega=w,dispersion_relative_defect=None if q==0. else abs(response(q/w)/(q*q)-1.)))
    openings=[]
    for c in (-1.,0.,1.):
        def qn(q):
            V=np.hypot(q-Q*c,Q*np.sqrt(1-c*c))
            return pole(q)-E-V
        z=brentq(qn,0.,.999*kL,xtol=2e-16)
        if qn(.999*z)<=0 or qn(1.001*z)>=0:
            raise ArithmeticError('Longitudinal past-null opening did not bracket')
        openings.append(dict(cosine=c,k=z,omega=pole(z)))
    return dict(grid=g,omega_plasma=wp,timelike_longitudinal_endpoint_k=kL,
                same_gap_small_k_coefficient=slope,small_k_relative_difference=abs(accepted/slope-1.),
                stable_longitudinal_branch_radial_relative_difference=cuterr,
                past_null_absorption_openings=openings,rows=rows)


def run():
    bag=[]
    p=Path(__file__).parent/'results'
    for g in ('fine','third'):
        with gzip.open(p/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:d=json.load(f)
        bag.append(evaluate(g,d))
    return dict(scope='Collisionless isotropic longitudinal pole only, fixed charged spectra and mirrored anti proxy; no pole emission rate.',data=bag)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
