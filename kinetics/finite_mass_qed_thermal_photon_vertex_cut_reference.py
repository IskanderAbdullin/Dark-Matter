"""The symmetric bare-photon thermal cut of the unequal-mass scalar vertex."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from early_pole_decay_mode_clock_reference import Q, intervals
from early_pole_hard_yukawa_photoproduction_reference import MET, dot, slash
from finite_mass_qed_soft_subtraction_reference import angular_lab


def projected(P,b,k,M,m):
    d=dot(P,b)
    dP=dot(P,k)
    db=dot(b,k)
    k2=dot(k,k)
    spin=(M*db+m*dP)/(d+M*m)
    return 4*(d-dP-db+k2)+2*(M+m)*spin


def trace(P,b,k,M,m):
    I=np.eye(4,dtype=complex)
    X=slash(b-k)+m*I
    Y=slash(P-k)+M*I
    G=slash(np.eye(4))
    V=sum(x*g@X@Y@g for x,g in zip(MET,G,strict=True))
    B=slash(b)+m*I
    A=slash(P)+M*I
    return (np.trace(B@V@A)/np.trace(B@A)).real


def samples(P,b,M,m,v):
    out=[]
    for direction in (np.array([1.,0.,0.]),np.array([0.,0.,1.]),
                      np.array([0.,0.,-1.]),np.r_[.6,0.,.8],np.r_[0.,.8,-.6]):
        D=np.r_[1.,v*direction]
        E=dot(P,D)
        E2=dot(b,D)
        d=dot(P,b)
        soft=2*d/(E*E2)
        # Compute k² in terms of v to keep the exactly null branch null.
        e=[]
        for rho in (1e-4,1e-5,1e-6):
            zz=1-v*v
            N=[projected(P,b,rho*D*sig,M,m) for sig in (-1.,1.)]
            exact=[trace(P,b,rho*D*sig,M,m) for sig in (-1.,1.)]
            vals=[N[j]/((rho*zz-2*sig*E)*(rho*zz-2*sig*E2))
                  for j,sig in enumerate((-1.,1.))]
            e.append(dict(plasma_photon_energy=rho,relative_to_eikonal=float(sum(vals)/soft-1),
                          max_direct_spin_trace_defect=float(max(abs(np.asarray(N)/exact-1)))))
        out.append(dict(direction=direction.tolist(),cases=e))
    return out


def thermal_angle_integral(P,b,M,m,rho,v):
    d=dot(P,b)
    Ep,Pz=P[[0,3]]
    Eb,Bx,Bz=b[[0,1,3]]
    # At each polar photon angle average the full shifted denominator
    # numerically over azimuth. For the null test use exact k²=0.
    phi=np.cos((np.arange(72)+.5)*2*np.pi/72)
    def inner(x):
        z=1-x
        Ep0=M*M/(Ep+Pz)+Pz*x
        ep=Ep0+(1-v)*(Pz*z)
        eb=(Eb-Bz)+Bz*x+(1-v)*Bz*z-v*Bx*np.sqrt(x*(2-x))*phi
        k2=rho*rho*(1-v*v)
        spin=(M*eb+m*ep)/(d+M*m)
        Np=4*(d-rho*(ep+eb)+k2)+2*(M+m)*rho*spin
        Nm=4*(d+rho*(ep+eb)+k2)-2*(M+m)*rho*spin
        a=rho*(1-v*v)
        full=Np/((a-2*ep)*(a-2*eb))+Nm/((a+2*ep)*(a+2*eb))
        leading=2*d/(ep*eb)
        return .5*np.mean(full-leading)
    x=np.r_[0.,np.geomspace(1e-8,.02,12),2.]
    return float(sum(quad(inner,a,b,epsabs=3e-13,epsrel=1e-9)[0]
                     for a,b in zip(x[:-1],x[1:],strict=True)))


def run():
    out=[]
    for g in ('fine','third'):
        p=Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
        with gzip.open(p) as f:d=json.load(f)
        M=d['x_end'];m=M*d['mass_ratio'];s=d['checkpoints'][-1]['scalar_pole_over_T']
        Es,lo,hi,*_=intervals(M,s,m,Q)
        Ep=(lo+hi)/2
        Pz=np.sqrt((Ep-M)*(Ep+M))
        c=(Ep*Es-(M*M+s*s-m*m)/2)/(Pz*Q)
        P=np.array([Ep,0.,0.,Pz])
        b=P-np.r_[Es,Q*np.sqrt((1-c)*(1+c)),0.,Q*c]
        if abs(dot(b,b)-m*m)>2e-15:
            raise ValueError('Archive recoil controls are not on shell.')
        massless=samples(P,b,M,m,1.)
        timelike=samples(P,b,M,m,.9)
        allrows=[a for s in (massless,timelike) for b0 in s for a in b0['cases']]
        if (max(abs(a['relative_to_eikonal']) for x in massless for a in x['cases'])>1e-13
                or max(a['max_direct_spin_trace_defect'] for a in allrows)>2e-9):
            raise ValueError('Projected cut or bare photon limit failed.')
        no=[thermal_angle_integral(P,b,M,m,r,1.) for r in (1e-4,1e-5)]
        stress=[thermal_angle_integral(P,b,M,m,r,.9) for r in (1e-4,1e-5)]
        ratio=stress[0]/stress[1]
        if max(abs(x) for x in no)>5e-13 or abs(ratio/100-1)>3e-4:
            raise ValueError('Off-null stress did not distinguish the exactly null identity.')
        v0,_,_=angular_lab(M,m,s,Q)
        out.append(dict(grid=g,bare_photon_integrated_leading=float(v0+2),null_point_tests=massless,
                        null_integrated_differences=no,timelike_toy_point_tests=timelike,
                        timelike_integrated_differences=stress,timelike_ratio=ratio))
    return dict(scope='Photon cut of scalar vertex only, Feynman gauge and vacuum external spinors. No finite thermal charged cut, residues or rate.',
                rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
