"""Entropy-checked longitudinal moment equations for a static thermal control."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fd_number_heat_state_reference import evaluate
from positive_fd_portal_quadrature_reference import PositivePortal


def model(T,xi,portal,receipt,k,proxy=None):
    h=evaluate(T,0.)
    n,d,chi,Cv,n_eta=(h[q] for q in
                      ("n","rho","chemical_energy","cv_fixed_number","dn_deta"))
    d+=h["pressure"]
    v=xi*d
    Cv_v=3*v/T
    _,_,_,m0,m1,m2,_=portal.moments(T,T,0.)
    Kq=receipt["timelike_heat"]+receipt["spacelike_heat"]
    Kp=receipt["timelike_bulk"]+receipt["spacelike_bulk"]
    gamma=Kq/d
    if proxy=="heat_as_momentum":
        Kp=Kq
    if proxy=="ignore_elastic_momentum":
        Kp=receipt["timelike_bulk"]

    # Coordinates: relative n_D, relative enthalpy-scaled rho_D and
    # rho_visible, followed by the longitudinal peculiar velocities.
    basis=np.eye(5)
    dn,dd,dv,ud,uv=(s*z for s,z in zip((n,d,v,1.,1.),basis,strict=True))
    td=(dd-chi*dn)/Cv
    tv=dv/Cv_v
    eta=dn/n_eta-chi*td/T**2
    pd=(d/T)*td+n*T*eta
    pv=dv/3
    number=2*m0*eta+(m1/T)*(td-tv)/T
    heat=2*m1*eta+Kq*(td-tv)/T
    drag=Kp*(ud-uv)
    L=np.array([-(1j*k*n*ud+2*number/gamma)/n,
                -(1j*k*d*ud+heat/gamma)/d,
                (-1j*k*v*uv+heat/gamma)/v,
                -(1j*k*pd+drag/gamma)/d,
                (-1j*k*pv+drag/gamma)/v])
    H=(np.outer(dn,dn)/n_eta + Cv*np.outer(td,td)/T**2
       + Cv_v*np.outer(tv,tv)/T**2 + (d*np.outer(ud,ud)+v*np.outer(uv,uv))/T)
    a,b=2*eta,(td-tv)/T**2
    D=(m0*np.outer(a,a)+m1*(np.outer(a,b)+np.outer(b,a))
       +(m2+T*receipt["spacelike_heat"])*np.outer(b,b)
       +(Kp/T)*np.outer(ud-uv,ud-uv))/gamma
    metric=H@L+L.conj().T@H+2*D
    closure=np.max(abs(metric))/(np.max(abs(H))*(1+k))
    if closure>3e-11 or np.linalg.eigvalsh(D).min()<-1e-12*np.max(abs(D)):
        raise ValueError("Streaming or exchange failed the independent entropy budget.")
    eig=np.linalg.eigvals(L)
    if eig.real.max()>3e-9*max(1.,k):
        raise ValueError("Unexpected growing local perturbation.")
    lock=(d+v)/(T*(Cv+chi**2*n_eta/T**2+Cv_v))
    temperature_in_free_wave=(d-chi*n)/Cv
    chemically_frozen_speed_squared=((d/T)*temperature_in_free_wave+
          n*T*(n/n_eta-chi*temperature_in_free_wave/T**2))/d
    return L,H,eig,closure,float(Kp/Kq*(1+1/xi)),float(lock),float(chemically_frozen_speed_squared)


def check_thermodynamics(T):
    m=evaluate(T,0.)
    h=m['rho']+m['pressure']
    n,c,chi,b=(m[k] for k in ('n','cv_fixed_number','chemical_energy','dn_deta'))
    theta,z=.31,-.21
    ndiff=b*(z+chi*theta/T)
    target=np.array([ndiff,chi*ndiff+c*T*theta,h*theta+n*T*z])
    e=2e-5
    vals=[]
    for g in (-2.,-1.,1.,2.):
        q=evaluate(T*(1+g*e*theta),g*e*z,nodes=174)
        vals.append(np.array([q[k] for k in ('n','rho','pressure')]))
    finite=(vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*e)
    err=float(np.max(abs(finite/target-1)))
    if err>1e-7:
        raise ValueError("Differential thermodynamic closure did not match finite EOS probes.")
    return err


def run(receipt):
    portal=PositivePortal((13,36,26),tail=65.)
    out=[]
    for source in receipt["rows"]:
        T=source["T"]
        thermo=check_thermodynamics(T)
        for xi in (1.,10.):
            rec=[]
            for scenario in (None,"heat_as_momentum","ignore_elastic_momentum"):
                L,H,eig,clos,relative,lock,free=model(T,xi,portal,source,0.,scenario)
                dT=np.array([0.,0.,0.,1.,1.])
                if np.max(abs(L@dT))>2e-13:
                    raise ValueError("Uniform common boost did not conserve momentum.")
                for conserved in (np.array([0.,1.,xi,0.,0.]),np.array([0.,0.,0.,1.,xi])):
                    if np.max(abs(conserved@L))>2e-12*xi:
                        raise ValueError("Uniform energy or momentum failed to close.")
                if min(abs(eig+relative))>1e-10:
                    raise ValueError("Exact relative drift damping did not replay.")
                # Both fluids start with opposite momentum, no density perturbation.
                start=np.array([0.,0.,0.,1.,-1./xi])
                measurements=[]
                for k in (.003,.03,.3,1.,3.,10.):
                    a,h,ev,ce,*_=model(T,xi,portal,source,k,scenario)
                    init=float(np.real(start@h@start))
                    previous=1.
                    for step in (.01,.1,1.,3.,10.):
                        end=expm(a*step)@start
                        now=float(np.real(end.conj()@h@end)/init)
                        if not 0.<=now<=previous+2e-9:
                            raise ValueError("Entropy norm increased in local time evolution.")
                        previous=now
                    end=expm(a)@start
                    measurements.append(dict(k_over_heat_rate=k,ep=[(float(z.real),float(z.imag)) for z in ev],
                                             relative_velocity_at_one_heat_time=abs(end[3]-end[4])/abs(start[3]-start[4]),
                                             metric_closure=ce))
                k=.00002
                u=model(T,xi,portal,source,k,scenario)[2]
                speed=max(abs(u.imag))/k
                if abs(speed/np.sqrt(lock)-1)>1e-3:
                    raise ValueError("Low-frequency locked sound speed did not match the LTE state.")
                k_fast=2e6
                high=model(T,xi,portal,source,k_fast,scenario)[2]
                cfast=np.sort(np.abs(high.imag)/k_fast)[[1,3]]
                if np.max(abs(cfast-np.sort(np.sqrt([free,1/3]))))>1e-4:
                    raise ValueError("High-frequency waves did not approach separate frozen sound speeds.")
                rec.append(dict(scenario=scenario or "full_momentum",uniform_relative_velocity_decay_per_heat_rate=relative,
                                locked_sound_speed_squared=lock,low_k_speed=float(speed),finite_wavelength_checks=measurements,
                                uncoupled_dark_sound_speed_squared=free,high_k_speeds=cfast.tolist()))
            out.append(dict(T=T,visible_to_dark_enthalpy_ratio=xi,finite_thermodynamic_derivative_error=thermo,controls=rec))
    return dict(scope="Static linear thermal two-fluid model; no gravity/FRW/scalar-transfer prediction.", rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt')
    a=p.parse_args()
    print(json.dumps(run(json.loads(Path(a.receipt).read_text())),indent=2,allow_nan=False))
