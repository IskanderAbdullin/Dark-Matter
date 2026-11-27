"""A conditional photon-temperature zero with the charged/scalar archive fixed."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_decay_mode_clock_reference import Q
from finite_mass_qed_crossed_annihilation_source_reference import inverse_tail_bound, value


def run_level(data,n,na):
    measure=[]
    frozen=value(data,1.,68.,n,-1,na,4,capture=measure)
    q,A,B,_=np.concatenate(measure,axis=0).T
    fs=frozen['scalar_f']
    pp=max(data['endpoint_parent_occupations'])
    pb=max(data['endpoint_fermion_occupations'])
    pointwise_sign_margin=1-pp-pb-pp*pb/fs
    if pointwise_sign_margin<=0. or min(A.min(),B.min())<0.:
        raise ArithmeticError('Pointwise thermal derivative sign condition failed')
    delta=(1+1/fs)*A-B
    if delta.max()>5e-16:
        raise ArithmeticError('Quadrature breached the pointwise sign bound')
    vacuum=(1+1/fs)*A.sum()
    def evaluate(T):
        if T==0.:
            return float(vacuum)
        n0=np.exp(-q/T)/(-np.expm1(-q/T))
        return float(vacuum+n0@delta)
    root=brentq(evaluate,0.,1.,xtol=2e-13)
    if abs(evaluate(1.)-frozen['net']/fs)>2e-12:
        raise ArithmeticError('Reweighted measure differs from independently summed collision clock')
    n0=np.exp(-q/root)/(-np.expm1(-q/root))
    grad=float((q*n0*(1+n0))@delta/root**2)
    if grad>=0.:
        raise ArithmeticError('Incorrect analytic T derivative at the root')
    h=1e-5
    check=abs((evaluate(root+h)-evaluate(root-h))/(2*h*grad)-1.)
    if check>5e-8:
        raise ArithmeticError(('T derivative quadrature failed',check))
    n1=np.exp(-q)/(-np.expm1(-q))
    loss_1=B*n1
    delta_1=-delta*n1
    c=[]
    for v in (2.,3.,5.,8.,16.,32.):
        mask=q<v
        c.append(dict(k_below=v,
                      fraction_of_T1_inverse_loss=float(loss_1[mask].sum()/loss_1.sum()),
                      fraction_of_T1_net_thermal_decrease=float(delta_1[mask].sum()/delta_1.sum())))
    return dict(T_zero=float(root),T_derivative=float(grad),T_derivative_relative_check=float(check),
                pointwise_negative_margin=float(pointwise_sign_margin),
                T0_net_per_fs=float(vacuum),T1_net_per_fs=evaluate(1.),
                equilibrium_temperature_is_not_being_inferred=True,
                tail_bound_on_net_per_fs_at_root=inverse_tail_bound(data,root,68.),
                rows=c)


def run():
    rows=[]
    for g in ('fine','third'):
        with gzip.open(Path(__file__).parent/'results'/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as f:
            d=json.load(f)
        low,high=(run_level(d,*mesh) for mesh in ((3,4),(6,8)))
        if abs(high['T_zero']/low['T_zero']-1)>2e-4:
            raise ArithmeticError('Frozen temperature zero did not converge')
        rows.append(dict(grid=g,low=low,high=high))
    return dict(scope='Fixed archived masses and charged/scalar histograms, equal anti-daughter proxy; only assumed photon Bose temperature varies; one crossed contribution at fixed q.',q=Q,rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
