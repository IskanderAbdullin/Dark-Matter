"""Non-null charged absorption trace, physical support, and sequential-pole diagnostic."""
import gzip,json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_hard_yukawa_photoproduction_reference import adj,dot,slash,I,norm_spinors
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q,case


def setup(d,k,w,c):
    M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
    E=np.hypot(s,Q)
    q=np.array([E,Q*np.sqrt(1-c*c),0.,Q*c])
    K=np.array([w,0.,0.,k])
    Qp=q-K
    V=np.linalg.norm(Qp[1:])
    def fun(e,phi):
        axis=Qp[1:]/V
        cross=np.array([axis[2],0.,-axis[0]])
        pp=(e*Qp[0]-(M*M+dot(Qp,Qp)-m*m)/2)/V
        rad=(e-M)*(e+M)-pp*pp
        if rad<0:raise ValueError('No parent angle')
        P=np.r_[e,pp*axis+np.sqrt(rad)*(np.cos(phi)*cross+np.sin(phi)*np.array([0.,1.,0.]))]
        b=P-Qp
        return P,b
    return M,m,q,K,Qp,fun


def spin(P,b,K,M,m,eps):
    Str,B,T,L=map(slash,(P,b,K,eps))
    s=dot(K,K)
    A=s+2*dot(P,K)
    D=s-2*dot(b,K)
    V=(Str+T+M*I)@L/A+L@(B-T+m*I)/D
    return float(norm_spinors(P,b,M,m,V))


def run_one(d,g):
    rows=[]
    for k,w,c in [(.0008,.0002,0.),(.00128,.00064,-.5),(.00256,.00064,.5),(.00256,
            next(z['omega'] for z in case(d)['poles'] if z['k']==.00256),0.)]:
        M,m,q,K,T,geo=setup(d,k,w,c)
        iscut=dot(K,K)<0
        grid=np.linspace(.12,5.,620)
        phase_cases=[]
        ward=gauge=denom=0.
        for ph in (0.,1.07,2.4,3.1):
            def den(e):
                P,b=geo(e,ph)
                return np.array([dot(K,K)+2*dot(P,K),dot(K,K)-2*dot(b,K)])
            permitted=[]
            for e in grid:
                try:
                    P,b=geo(e,ph)
                except ValueError:
                    continue
                if b[0]>m:permitted.append(e)
            roots=[]
            for arm in (0,1):
                for xx,yy in zip(permitted[:-1],permitted[1:],strict=True):
                    if den(xx)[arm]*den(yy)[arm]<0:
                        z=brentq(lambda e:den(e)[arm],xx,yy,xtol=2e-14)
                        roots.append(dict(arm=arm,parent_E=z,denominators=den(z).tolist()))
            for e in permitted[::19]:
                P,b=geo(e,ph);de=den(e)
                denom=max(denom,abs(de.sum()-2*dot(q,K)),abs(dot(b,b)-m*m))
                if min(abs(de))>1.e-8:
                    V0=np.array([0.,1.,0.,0.])
                    V1=np.array([0.,0.,1.,0.])
                    t=spin(P,b,K,M,m,V0)+spin(P,b,K,M,m,V1)
                    if t<=0:raise ArithmeticError('Nonpositive transverse trace')
                    z=spin(P,b,K,M,m,V0+.8*K)+spin(P,b,K,M,m,V1-.4*K)
                    gauge=max(gauge,abs(z/t-1))
                    ward=max(ward,abs(spin(P,b,K,M,m,K)/t))
            phase_cases.append(dict(azimuth=ph,found_near_born_intermediate_on_shell_roots=roots))
        if gauge>3e-8 or ward>3e-8 or denom>1e-10:
            raise ArithmeticError((g,k,gauge,ward,denom))
        rows.append(dict(k=k,omega=w,vacuum_photon_invariant=dot(K,K),charged_transfer_invariant=dot(T,T),
                         analytic_sum_of_two_denominators=2*dot(q,K),
                         transverse_shift_gauge_relative_defect=gauge,ward_relative_defect=ward,
                         kinematics_and_invariant_denominator_defect=denom,azimuths=phase_cases))
    return dict(grid=g,checks=rows)


def run():
    out=[]
    for g in ('fine','third'):
        p=Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:d=json.load(f)
        out.append(run_one(d,g))
    return dict(scope='Bare spinor non-null-photon amplitude and location of kinematic bare denominators; no assignment of spectral rates, widths or subtraction.',data=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
