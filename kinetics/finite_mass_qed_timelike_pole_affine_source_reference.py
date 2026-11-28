"""Fixed-k, fixed-angle absorption gain/loss with a free pole occupation variable."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I, norm_spinors, slash
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q, Species
from finite_mass_qed_longitudinal_static_spectral_reference import DensityResponse


def vertex(P,B,K,m0,m1,e,pieces=False):
    """Use stable invariant denominators with equivalent bare numerators."""
    mu2=(K[0]-K[3])*(K[0]+K[3])
    vals=[]
    for p,m in ((P,m0),(B,m1)):
        q=np.linalg.norm(p[1:])
        eq=np.hypot(m,q)
        longitudinal=p[3]
        turn=p[1:3]@p[1:3]/(q+longitudinal) if longitudinal>0 else q-longitudinal
        stable=(m*m*K[0]*K[0]+mu2*q*q)/(eq*K[0]+q*K[3])+K[3]*turn
        direct=p[0]*K[0]-p[3]*K[3]
        if abs(stable-direct)>3e-10*stable:
            raise ArithmeticError(('Stable on-shell denominator mismatch',stable,direct,p,K))
        vals.append((mu2+2*stable,mu2-2*stable))
    Y=slash(e)
    left=(slash(P+K)+m0*I)@Y/vals[0][0]
    right=Y@(slash(B-K)+m1*I)/vals[1][1]
    return (left,right) if pieces else left+right


def setup(m0,m1,k,w,es,c):
    Q0=np.array([es,Q*np.sqrt(1-c*c),0.,Q*c])
    K=np.array([w,0.,0.,k])
    T=Q0-K
    V=np.linalg.norm(T[1:])
    t=T[0];s=(t-V)*(t+V);A=m0*m0+s-m1*m1
    if t>=0. or s>=0.:
        raise ValueError('Pilot expects above the past-null, space-like transfer opening')
    Emin=(-A*t+V*np.sqrt(A*A-4*m0*m0*s))/(-2*s)
    axis=T[1:]/V
    cross=np.array([axis[2],0.,-axis[0]])
    def shell(E,angle):
        par=(E*t-A/2)/V
        r2=(E-m0)*(E+m0)-par*par
        if r2<0.:raise ArithmeticError('Outside exact shell during quadrature')
        P=np.r_[E,par*axis+np.sqrt(r2)*(np.cos(angle)*cross+np.sin(angle)*np.array([0.,1.,0.]))]
        return P,P-T
    return K,t,V,Emin,shell


def run_one(d,S,g,name,k,c,n,az):
    m0,m1=S[0].m,S[1].m
    es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    if name=='T':
        def inv(w):
            return (w-k)*(w+k)-sum(x.field(k/w)[0] for x in S)
        w=brentq(inv,np.nextafter(k,np.inf),k+.006,xtol=5e-17)
        response_der=sum(x.field(k/w)[1] for x in S)
        residue=1/(2*w+k/(w*w)*response_der)
        pol=np.eye(4)[1:3]
    else:
        w,residue=DensityResponse(S).pole(k)
        model=DensityResponse(S)
        inv=lambda a:k*k-model.retarded(a/k)[0]
        pol=np.eye(4)[:1]
    h=(w-k)*4e-5
    central=((inv(w+h)-inv(w-h))/(2*h))*residue
    derivative_defect=abs(central-1)
    if derivative_defect>5e-8:
        raise ArithmeticError(('Physical residue derivative',g,name,k,derivative_defect))
    K,t,V,lo,shell=setup(m0,m1,k,w,es,c)
    hi=max(S[0].E[-1],S[1].E[-1]+t)
    if not lo<hi:raise ArithmeticError(('No archived support',g,name,k,lo,hi))
    ed=np.unique(np.r_[lo,hi,S[0].E,S[1].E+t])
    ed=ed[(ed>=lo)&(ed<=hi)]
    fext=[np.r_[0.,x.f,0.] for x in S]
    x,wei=roots_legendre(n)
    stats=np.zeros(3); thermal=np.zeros(2); ward=0.; ward_weight=0.; trace_weight=0.; single_ward=0.; gauge_weight=0.; gauge_ref=0.
    for a,b in zip(ed[:-1],ed[1:],strict=True):
        mid=(a+b)/2
        p=fext[0][np.searchsorted(S[0].E,mid)]
        y=fext[1][np.searchsorted(S[1].E,mid-t)]
        if p==y==0.:
            continue
        for z,u in zip(mid+(b-a)*x/2,wei*(b-a)/2,strict=True):
            trace=[]
            for angle in (np.arange(az)+.5)*2*np.pi/az:
                P,B=shell(z,angle)
                den=sum(norm_spinors(P,B,m0,m1,vertex(P,B,K,m0,m1,e)) for e in pol)
                if den<=0:raise ArithmeticError('Nonpositive squared radiative matrix')
                trace.append(den)
            P,B=shell(z,0.198)
            a1,a2=vertex(P,B,K,m0,m1,K/w,pieces=True)
            born=norm_spinors(P,B,m0,m1,I/w)
            single_ward=max(single_ward,*(norm_spinors(P,B,m0,m1,X)/born for X in (a1-I/w,a2+I/w)))
            W=norm_spinors(P,B,m0,m1,a1+a2)
            ward=max(ward,W/np.mean(trace))
            if name=='L':
                density=norm_spinors(P,B,m0,m1,vertex(P,B,K,m0,m1,np.eye(4)[0]))
                s=np.array([0.,0.,0.,-k/w])
                j=norm_spinors(P,B,m0,m1,vertex(P,B,K,m0,m1,s))
                gauge_weight+=u*abs(density-j)
                gauge_ref+=u*(density+j)/2
            ward_weight+=u*W
            trace_weight+=u*np.mean(trace)
            W=u*np.mean(trace)*residue/V
            stats+=W*np.array([p*(1-y),y*(1-p),p-y])
            P0,B0=(expit((.13-E)/.8) for E in (z,z-t))
            thermal+=W*np.array([P0*(1-B0),B0*(1-P0)])
    if min(stats[:2])<=0 or ward>2e-11 or ward_weight/trace_weight>5e-14 or single_ward>2e-14:
        raise ArithmeticError(('Gain/loss or Ward gate',g,name,k,c,stats,ward,ward_weight/trace_weight,single_ward))
    if abs(stats[0]-stats[1]-stats[2])>5e-13*max(stats[:2]):
        raise ArithmeticError('Independent signed Pauli integration failed')
    nS=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    slope=stats[2]+stats[0]/nS
    ncrit=stats[1]/slope if slope>0. else None
    nthermal=1/np.expm1(es/.8)
    nph=thermal[1]/((1+1/nthermal)*thermal[0]-thermal[1])
    thermal_occupation_error=abs(nph*np.expm1(w/.8)-1)
    gauge=gauge_weight/gauge_ref if name=='L' else None
    if thermal_occupation_error>2e-11 or (gauge is not None and gauge>2e-7):
        raise ArithmeticError(('Thermal neutral or independent density projection gate',thermal_occupation_error,gauge))
    pnorm=np.array([stats[0],stats[1],stats[2]])
    return dict(grid=g,polarization=name,k=k,omega=w,c=c,mesh=[n,az],raw_kernel=pnorm.tolist(),
                support_parent_energy=lo,pole_residue_omega=residue,pole_residue_derivative_mesh_error=derivative_defect,
                pole_population_neutral=ncrit,largest_ward_squared=float(ward),weighted_ward_squared=float(ward_weight/trace_weight),single_diagram_ward_squared=float(single_ward),
                common_thermal_relative_defect=abs(thermal[1]/thermal[0]/np.exp(t/.8)-1),
                thermal_neutral_I2_relative_defect=thermal_occupation_error,longitudinal_density_vs_current_gauge_defect=gauge)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--refine',action='store_true')
    o=p.parse_args()
    path=Path(__file__).parent/'results'
    with gzip.open(path/f'early-pole-spectrum-{o.grid}-m06-2026-11-26.json.gz') as fd:
        d=json.load(fd)
    M=d['x_end']
    S=[Species(m,d[f'endpoint_{j}_momentum_edges'],d[f'endpoint_{j}_occupations'])
       for m,j in ((M,'parent'),(M*d['mass_ratio'],'fermion'))]
    for name,k,c in (('T',.0015,0.),('T',.00256,0.),('L',.0012,0.),('L',.00256,0.)):
        r=run_one(d,S,o.grid,name,k,c,*( (6,12) if o.refine else (3,6)))
        print(json.dumps(r,allow_nan=False),flush=True)
