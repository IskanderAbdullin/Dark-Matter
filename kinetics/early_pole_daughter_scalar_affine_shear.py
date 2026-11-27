"""Instantaneous shear of a frozen isotropic scalar state, not a persisted response."""
import argparse
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from early_pole_daughter_bare_scalar_stress_control import state
from early_pole_daughter_retained_kinetic_stress import build, r


def deformed(q, n, lam, bare, eta, angles):
    mu, weights = leggauss(angles)
    q2 = q[:, None]**2
    lateral = np.exp(-eta)*q2*(1-mu**2)
    axial = np.exp(2*eta)*q2*mu**2
    nangle = n[:,None] * .5*weights
    p2 = lateral+axial
    summed = lambda x: float(np.sum(nangle*x))
    eq = lambda m: m*m-bare*bare-.5*lam*summed(1/np.sqrt(p2+m*m))
    hi = 2*(bare+(.5*lam*float(sum(n)))**(1/3))
    m = brentq(eq, max(1e-14,bare),hi,xtol=8e-18)
    E = np.sqrt(p2 + m*m)
    V = (m*m-bare*bare)**2/(2*lam)
    return np.array([m, summed(E)-V, summed(p2/(3*E))+V,
                     summed((axial-.5*lateral)/E)])


def check(q, n, lam, bare):
    z=state(q,n,lam,bare)
    m=z['effective_gap']
    E=np.hypot(q,m)
    pkin=float(n@(q*q/(3*E)))
    response=float(n@(q*q/E*(1-(q/E)**2/5)))
    if not 2.4*pkin-8*np.spacing(pkin)<=response<=3*pkin+8*np.spacing(pkin):
        raise ValueError('Analytic affine-shear bounds did not close.')
    # Self-consistent anisotropy finite differences use a separate gap solver
    # for each shape. Pressure is per principal axis, transverse axes equal.
    checks=[]
    for order in (8,16):
        h=.002
        f=np.array([deformed(q,n,lam,bare,j*h,order) for j in (-2,-1,0,1,2)])
        d=(f[0]-8*f[1]+8*f[3]-f[4])/(12*h)
        d2=(-f[0]+16*f[1]-30*f[2]+16*f[3]-f[4])/(12*h*h)
        if (abs(d[0])>3e-9*m or abs(d[1])>3e-9*z['energy'] or abs(d[2])>3e-9*z['pressure']
                or abs(f[2][2]-z['pressure'])>2e-12*z['pressure']):
            raise ValueError('Scalar singlets did not have vanishing first-order shear.')
        if abs(d[-1]-response)>4e-9*max(pkin,1e-25):
            raise ValueError('Anisotropic stress tensor did not agree with its analytic slope.')
        if abs(d2[1]-response)>2e-7*max(pkin,1e-12*z['energy'])+1e-9*z['energy']:
            raise ValueError('Work curvature did not agree with stress derivative.')
        checks.append(dict(order=order, response=float(d[-1]), energy_curvature=float(d2[1]),
                           first_energy_derivative=float(d[1]), first_gap_derivative=float(d[0])))
    return dict(w=z['pressure_to_energy'], bare_mass=bare, effective_gap=m, affine_shear=float(response),
                affine_shear_over_pressure=response/z['pressure'],
                affine_shear_over_energy=response/z['energy'], kinetic_pressure_fraction=pkin/z['pressure'],
                finite_checks=checks)


def paired_response(q, left, right, a, b, cold_modes):
    # Split a paired difference without subtracting near-identical stress totals.
    m,H=a['effective_gap'],b['effective_gap']
    A,B=np.hypot(q,m),np.hypot(q,H)
    u,v=1/A,1/B
    h=q*q
    inv_gap=-(H-m)*(H+m)/(A*B*(A+B))
    response_gap=h*inv_gap*(1 - h*(u*u+u*v+v*v)/5)
    response_count=h*u*(1-h*u*u/5)
    count=float((right-left)@response_count)
    pole_cold=float(right[:cold_modes]@response_gap[:cold_modes])
    pole_spect=float(right[cold_modes:]@response_gap[cold_modes:])
    native_cold=float(left[:cold_modes]@response_count[:cold_modes])
    net=count+pole_cold+pole_spect
    if abs(net-(b['affine_shear']-a['affine_shear']))>32*np.spacing(a['affine_shear']):
        raise ValueError('Subtracted shear did not reproduce.')
    return dict(absolute_paired_response=net, relative_paired_response=net/a['affine_shear'],
                native_retained_response_fraction=native_cold/a['affine_shear'],
                retained_count_shape_term=count, retained_pole_term=pole_cold, exterior_pole_term=pole_spect)


def run(archive, arms, receipt, angular=None):
    C,hard=build(archive,angular)
    if receipt['hard_depth']!=hard or receipt['cold_modes']!=len(C.q):
        raise ValueError('Unmatched scalar arm geometry.')
    q=np.r_[C.q,C.spq]
    out=[]
    nums=[]
    for t,d,tag in zip(arms,receipt['arms'],('native','local-daughter'),strict=True):
        if t.get('prepared_charged_swap','native')!=tag or d['arm']!=tag or t.get('aperture_subdivision',1)!=hard:
            raise ValueError('Unmatched terminal arm.')
        if abs(t['energy_defect'])>=3e-12 or abs(t['entropy_defect'])>=3e-11:
            raise ValueError('Original gates not met.')
        np.testing.assert_array_equal(t['cold_mode_count_history'][0],C.ns)
        np.testing.assert_array_equal(t['cold_momentum_nodes_over_T'],C.q)
        final=np.asarray(t['cold_mode_count_history'][-1])
        if not np.all(np.isfinite(final)&(final>0)):
            raise ValueError('Terminal populations invalid.')
        n=np.r_[final,C.spect]
        nums.append(n)
        masses=[('massless',0.)]+[(k,d['cases'][k]['bare_mass']) for k in d['cases']]
        known={'massless':d['baseline']['w'],**{k:m['w'] for k,m in d['cases'].items()}}
        z={k:check(q,n,C.lam,b) for k,b in masses}
        if any(abs(v['w']-known[k])>2e-12 for k,v in z.items()):
            raise ValueError('Pressure ratio did not match the homogeneous reference.')
        out.append(dict(arm=tag,cases=z))
    # Mass-matched pairs: each threshold receipt generally assigns slightly
    # different b to its two arms. Compare both at exactly the native b.
    paired={}
    for k,a in out[0]['cases'].items():
        b=check(q,nums[1],C.lam,a['bare_mass'])
        paired[k]=paired_response(q,*nums,a,b,len(C.q))
    zero=check(np.array([0.]),np.array([.75]),2.,.5)
    if zero['affine_shear']!=0.:
        raise ValueError('An exact zero-momentum sample acquired affine shear.')
    return dict(scope='Instantaneous volume-preserving affine momentum deformation; no dynamical anisotropic-stress closure.',
                hard_depth=hard,cold_modes=len(C.q),zero_momentum_fixture=zero,arms=out, paired_at_native_bare_mass=paired)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','native','local','homogeneous_receipt'):
        p.add_argument(k)
    p.add_argument('--angular-mesh')
    a=p.parse_args()
    print(json.dumps(run(a.archive,[r(a.native),r(a.local)],r(a.homogeneous_receipt),
                         r(a.angular_mesh) if a.angular_mesh else None),indent=2,allow_nan=False))
