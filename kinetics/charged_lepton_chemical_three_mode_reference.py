"""Passive 3-mode reduction of the local finite-chemistry acoustic control."""
import argparse
import json
from pathlib import Path

import numpy as np

from charged_lepton_portal_local_longitudinal_reference import model
from positive_fd_portal_quadrature_reference import PositivePortal


def isolated_peak(c0,c1):
    r=(c0/c1)**2
    if not 1/9<r<1:
        raise ValueError("Peak formula assumes a complex acoustic pair at all k>0.")
    q=r**.25
    sqrt=np.sqrt(r)
    d=(1-r)/(1+sqrt)
    expected=d/np.sqrt(4*r-d*d)
    def direct(k):
        roots=np.roots([1.,1.,k*k,r*k*k])
        acoustic=roots[int(np.argmax(roots.imag))]
        if acoustic.imag<=0. or acoustic.real>=0.:
            raise ValueError("Isolated cubic has no damped acoustic pair on this fixture.")
        return float(-acoustic.real/acoustic.imag)
    comparisons=[direct(v*q) for v in (.8,1.,1.25)]
    if abs(comparisons[1]/expected-1)>1e-8 or not max(comparisons[0],comparisons[2])<comparisons[1]:
        raise ValueError("Analytic peak disagrees with independent cubic roots.")
    real_root=-sqrt
    if abs(np.polyval([1.,1.,q*q,r*q*q],real_root))>3e-13:
        raise ValueError("The real chemical root does not close at the maximum.")
    return dict(squared_sound_speed_ratio=float(r),peak_k_over_gamma=float(q/c1),
                peak_attenuation_to_angular_frequency=float(expected),
                amplitude_loss_fraction_per_cycle=float(-np.expm1(-2*np.pi*expected)),
                direct_root_ratios_at_relative_k=[.8,1.,1.25],direct_attenuation_to_frequency=comparisons)


def one(portal,row,ref):
    T=row['T']
    ratio=ref['enthalpy_ratio']
    L,H,*_=model(T,ratio,portal,row,0.)
    E=model(T,ratio,portal,row,0.,elastic_only=True)[0]
    F=model(T,ratio,portal,row,1.,elastic_only=True)[0]
    w,V=np.linalg.eigh(H)
    W=(V*np.sqrt(w))@V.T
    J=(V/np.sqrt(w))@V.T
    trans=lambda X:W@X@J
    scale=1+row['timelike_heat']/row['spacelike_heat']
    A=-trans(E).real
    C=-scale*trans(L).real-A
    B=-trans(F-E).imag
    e,U=np.linalg.eigh((A+A.T)/2)
    S,R=U[:,:3],U[:,3:]
    c0,c1=(ref[key] for key in ('full_equilibrium_sound_speed','frozen_sound_speed'))
    alpha=ref['projected_chemical_rate_per_timelike_weight']
    residue=(c1**2-c0**2)/(2*alpha)
    if abs(residue/ref['projected_singular_attenuation_weight']-1)>1e-10:
        raise ValueError("Bare cubic did not reproduce the independently computed damping pole.")
    out=[]
    for case in ref['sweep']:
        eps=case['timelike_weight']
        D=A+eps*C
        dec,V=np.linalg.eigh((D+D.T)/2)
        if dec[2]<=0. or max(abs(dec[:2]))>2e-13:
            raise ValueError("Full loss did not retain the common energy and momentum nulls.")
        probe=[]
        for src in case['frequency_path']:
            k=src['wave_number']
            x=k/(eps*alpha)
            z=np.roots([1.,1.,c1*c1*x*x,c0*c0*x*x])
            a=z[int(np.argmin(z.imag))]
            cubic=(-a.imag/x,-a.real*eps*alpha/k**2)
            # Passive fast elimination: stationary fast response at fixed
            # slow coordinates. Keep the complete fast block including
            # its streaming. Any remainder is the omitted fast memory.
            M=-D-1j*k*B
            slave=np.linalg.solve(-(R.T@M@R),R.T@M@S)
            U0=S+R@slave
            eff=S.T@M@U0
            if np.max(abs(R.T@M@U0))>3e-13:
                raise ValueError("Fast coordinates did not satisfy their stationary equations.")
            eig,v=np.linalg.eig(eff)
            j=int(np.argmin(eig.imag))
            u=v[:,j]
            g=V[:,2:].T@(U0@u)
            loss=float(dec[2:]@abs(g)**2/np.vdot(u,u).real)
            if abs(loss+eig[j].real)>5e-13 or loss<0.:
                raise ValueError("Reduced model lost entropy passivity.")
            schur=(-eig[j].imag/k,loss/k**2)
            answer=src['phase_speed'],src['damping_over_k_squared']
            probe.append(dict(k_over_chem=src['k_over_slow_chemical_rate'],
                              cubic_phase_error_as_fraction_of_speed_gap=float(abs(cubic[0]-answer[0])/(c1-c0)),
                              cubic_damping_relative_error=float(abs(cubic[1]/answer[1]-1)),
                              schur_phase_error_as_fraction_of_speed_gap=float(abs(schur[0]-answer[0])/(c1-c0)),
                              schur_damping_relative_error=float(abs(schur[1]/answer[1]-1))))
        if max(r['schur_damping_relative_error'] for r in probe)>0.04:
            raise ValueError("Stationary elastic reduction exceeded the declared test domain.")
        out.append(dict(timelike_weight=eps,trace=probe))
    return dict(T=T,elastic_clock=1,analytic_cubic_damping_pole_residue=residue,
                isolated_chemical_peak=isolated_peak(c0,c1),sweep=out)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('portal_receipt')
    p.add_argument('crossover')
    a=p.parse_args()
    read=lambda s:json.loads(Path(s).read_text())['rows']
    isolated_peak(.35,.8)  # Independent positive fixture away from the near-conformal cases.
    original={r['T']:r for r in read(a.portal_receipt)}
    r=read(a.crossover)
    portal=PositivePortal((13,36,26),tail=65.)
    print(json.dumps(dict(scope="Three-mode static reduction bounded only in tested local domains.",
                          rows=[one(portal,original[v['T']],v) for v in r]),indent=2,allow_nan=False))
