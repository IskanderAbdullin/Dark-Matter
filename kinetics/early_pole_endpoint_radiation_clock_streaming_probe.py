"""Prescribed radiation-like conformal clock on frozen scalar distributions."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_frozen_ballistic_shape import eval_meas
from early_pole_endpoint_frozen_stress_sweep import nodal_frozen,weights
from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack


def average_speed(q,m,A):
    """Average q/sqrt(q^2+(m*a)^2) on 1<=a<=A; a is the prescribed clock."""
    if not np.isfinite(A) or not A>=1 or not np.isfinite(m) or m<0 or np.any(q<=0):
        raise ValueError('Positive support and a forward conformal interval required')
    if m==0:
        return np.ones_like(q)
    if A==1:
        return q/np.hypot(q,m)
    v=q*(np.arcsinh(m*A/q)-np.arcsinh(m/q))/(m*(A-1))
    if np.any(v<q/np.hypot(q,m*A)*(1-5e-14)) or np.any(v>q/np.hypot(q,m)*(1+5e-14)):
        raise ArithmeticError('Analytic path escaped monotone speed limits')
    return v


def clock_check(q,m,A):
    # Integrate in logarithmic a to keep the endpoints resolved for a large A.
    t,w=roots_legendre(64)
    x=np.exp(np.log(A)*(1+t)/2)
    expect=(q[:,None]/np.hypot(q[:,None],m*x))@(w*x*np.log(A)/(2*(A-1)))
    return float(np.max(abs(average_speed(q,m,A)/expect-1)))


def assemble(native,gauss):
    d,y=unpack(native,False,.05)
    M=d.state(.05,y)[0]
    probes=[('W1',[weights(d,y,n) for n in (96,192)])]
    for c in gauss:
        rep,f,k=nodal_frozen(c,d,M)
        w=[weights(d,f,n) for n in (96,192)]
        if abs(fsum(w[-1][1])/rep[0]-1)>3e-14:
            raise ArithmeticError('Count was not preserved in the projected spectrum')
        probes.append((f'G{k}',w))
    if [t for t,_ in probes]!=['W1','G3','G4']:
        raise ArithmeticError('Need the actual shared W1/G3/G4 baseline')
    return probes,M


def ratio_block(samples,u):
    H=np.array([fsum(p*v*v) for v,p in samples])
    b=u/np.sqrt(H)
    R=np.asarray([[eval_meas(v,p,t) for t in b] for v,p in samples])
    # Two symmetric paths in the two-dimensional table separate the
    # rate change at a fixed phase from the redefinition of phase itself.
    def split(i,j):
        shape=((R[j,i]-R[i,i])+(R[j,j]-R[i,j]))/2
        phase=((R[j,j]-R[j,i])+(R[i,j]-R[i,i]))/2
        D=R[j,j]-R[i,i]
        if abs(D-shape-phase)>3e-15:
            raise ArithmeticError('Symmetric decomposition does not close')
        return [float(D),float(shape),float(phase),float(R[j,0]-R[i,0])]
    return np.r_[R.ravel(),split(0,2),split(1,2)],b,H


def run(native,gauss):
    P,m=assemble(native,gauss)
    out=[]
    for A in (1.,2.,10.,100.,1000.):
        Z=[]
        errors=[]
        for _,rows in P:
            deck=[]
            for q,w,_ in rows:
                v=average_speed(q,m,A); p=w/fsum(w)
                deck.append((v,p))
            Z.append(deck)
            if A>1:
                x=rows[-1][0]
                errors.append(clock_check(x[::53],m,A))
        defect=max(errors,default=0.)
        if defect>5e-12:
            raise ArithmeticError('Analytic and independent flight integral disagree')
        rows=[]; reject=[]
        for u in (.003,.03,.3,1.,2.):
            a,b=[ratio_block([r[i] for r in Z],u) for i in (0,1)]
            err=float(np.max(abs(a[0]-b[0])))
            if err>max(3e-13,3e-9*np.max(b[0][:9])):
                reject.append(dict(u=u,absolute_quadrature_difference=err))
                continue
            rows.append(dict(u=u,raw_phase_for_each_rms=b[1].tolist(),second_moment_of_mean_speed=b[2].tolist(),
                             sampled_phase_matrix=b[0][:9].reshape(3,3).tolist(),
                             G4_minus_W1=b[0][9:13].tolist(),G4_minus_G3=b[0][13:].tolist(),
                             maximum_absolute_paired_defect=err))
        out.append(dict(A=A,analytic_vs_independent_clock_max_relative_defect=defect,accepted=rows,unresolved=reject))
    return dict(frozen_x=.05,initial_formal_mass=m,clock='uniform a in conformal time (prescribed radiation-like)',
                labels=['W1','G3','G4'],differences=['rms_matched','symmetric_shape','symmetric_phase','common_W1_raw_phase'],
                rows=out,scope='Independent free-characteristic with no metric forcing or kinetic coupling. Expansion ratio and phase are uncalibrated, not a cosmological transfer function.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--gauss',nargs=2,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(run(read(a.native),[read(f) for f in a.gauss]),indent=2,allow_nan=False))
