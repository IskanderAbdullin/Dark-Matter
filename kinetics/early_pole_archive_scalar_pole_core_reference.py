"""Signed local vacuum pole remainder under an explicitly chosen subtraction."""
import argparse
import gzip
import json

import numpy as np
from scipy.special import roots_legendre

from early_pole_archive_scalar_below_pole_finite_q_reference import integrate
from early_pole_native_aperture_witness_reference import GRID


def density(archive, q, order, x):
    row = integrate(archive,q,order,.9,probe_x=x,resolve_terms=True)
    return np.concatenate((np.array(row['term_gross_densities_per_reference_H']),
                           np.array(row['term_net_densities_per_reference_H'])[:,:,None]),axis=2)


def residual(radius, ratio, order, at_zero, evaluate, breaks=()):
    nodes,weights = roots_legendre(order)
    width = radius*ratio
    edges = np.unique([0.,radius,*[width*k for k in (1.,4.,16.) if 0.<width*k<radius],
                      *[b for b in breaks if 0.<b<radius]])
    h = np.diff(edges)[:,None]/2
    x = ((edges[1:]+edges[:-1])[:,None]/2+h*nodes).ravel()
    w = (h*weights).ravel()
    out = np.zeros((3,)+at_zero.shape[1:])
    for t,wt in zip(x,w,strict=True):
        plus,minus = evaluate(t),evaluate(-t)
        den = t*t+width*width
        even = (plus[0]-at_zero[0])+(minus[0]-at_zero[0])
        out += wt*np.array([even/den,t*(plus[1]-minus[1])/den,plus[2]+minus[2]])
    return out


def fixture():
    def function(x):
        return np.array([1+.2*x+3*x*x,4+.7*x,1+x])
    a = .083
    for r in (0.,.01,.1,1.):
        j = a-a*r*np.arctan(1/r) if r else a
        want = np.array([6*j,1.4*j,2*a])
        for n in (48,96):
            np.testing.assert_allclose(residual(a,r,n,function(0),function),want,rtol=3e-9)


def energy_integral(volume):
    p = np.cbrt(volume)
    m = .003
    e = np.hypot(p,m)
    return 3/8*(p*e*(2*p*p+m*m)-m**4*np.arcsinh(p/m))


def conserved_ramp(a, grid, fraction, subcells):
    p = np.asarray(a['endpoint_fermion_momentum_edges'])
    f = np.asarray(a['endpoint_fermion_occupations'])
    center = 1.4222222222222223 if grid=='fine' else 1.44
    i = int(np.argmin(abs(p-center)))
    assert abs(p[i]-center)<1e-10 and 0<i<len(p)-1
    v = p**3
    dx = fraction*min(v[i]-v[i-1],v[i+1]-v[i])
    left,right = v[i]-dx,v[i]+dx
    band = np.linspace(left,right,subcells+1)
    w = (np.arange(subcells)+.5)/subcells
    ramp = (1.-w)*f[i-1]+w*f[i]
    # Two shoulders preserve the original number and energy over the
    # union of the two neighboring native cells.
    ep = np.r_[v[:i],(v[i-1]+left)/2,band,(right+v[i+1])/2,v[i+1:]]
    out = np.r_[f[:i],f[i-1],ramp,f[i],f[i:]]
    dm,de = np.diff(ep),np.diff(energy_integral(ep))
    target = np.array([np.diff(v[i-1:i+2])@f[i-1:i+1],
                       np.diff(energy_integral(v[i-1:i+2]))@f[i-1:i+1]])
    area = slice(i-1,i+subcells+3)
    ids = [i-1,i+subcells+2]
    # The full-grid summation here also makes the earlier reference run
    # reproducible; independently check the affected two-cell union.
    all_moments = np.array([np.diff(v)@f,np.diff(energy_integral(v))@f])
    delta = np.linalg.solve(np.array([dm[ids],de[ids]]),all_moments-np.array([dm@out,de@out]))
    out[ids] += delta
    np.testing.assert_allclose([dm[area]@out[area],de[area]@out[area]],target,rtol=5e-12)
    np.testing.assert_allclose([dm@out,de@out],all_moments,rtol=5e-15)
    assert np.all((out>=0.)&(out<=1.))
    archive = dict(a,endpoint_fermion_momentum_edges=np.cbrt(ep).tolist(),
                   endpoint_fermion_occupations=out.tolist())
    controls = dict(ramp_momentum_bounds=np.cbrt([left,right]).tolist(),ramp_fraction=fraction,
                    ramp_subcells=subcells,occupation_adjustments_on_shoulders=delta.tolist())
    return archive,controls


def run(grid,q,radii,width_fractions,radial_orders,core_orders,ramp_fraction,subcells):
    with gzip.open(GRID/f'early-pole-spectrum-{grid}-m06-2026-11-26.json.gz','rt') as f:
        a = json.load(f)
    controls = {}
    if ramp_fraction:
        a,controls = conserved_ramp(a,grid,ramp_fraction,subcells)
    M,m = .05,.003
    ms = a['checkpoints'][-1]['scalar_pole_over_T']
    P = np.array(a['endpoint_fermion_momentum_edges'])
    E = np.hypot(P,m)
    x = np.abs((m*m+ms*ms+2*(E*np.hypot(q,ms)+np.array([-1.,1.])[:,None]*P*q))/(M*M)-1)
    breaks = np.unique(x[x<max(radii)])
    out = []
    for order in radial_orders:
        at_zero = density(a,q,order,0.)
        memo = {0.:at_zero}

        def evaluate(t):
            if t not in memo:
                memo[t] = density(a,q,order,t)
            return memo[t]

        for r in radii:
            for n in core_orders:
                values = []
                for g in width_fractions:
                    v = residual(r,g,n,at_zero,evaluate,breaks)
                    values.append(dict(gamma_over_R=g,term_residual=v.tolist(),combined_residual=v.sum(0).tolist(),
                                       shell_weight_fraction=(2/np.pi*np.arctan(1/g) if g else 1.)))
                out.append(dict(invariant_order=order,interval_radius=r,core_order=n,terms=values,
                                pole_coefficients=at_zero.tolist()))
    return dict(grid=grid,q=q,**controls,incoming_archive_topology_breaks=breaks.tolist(),results=out)


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    p.add_argument('--q',type=float,default=4.1e-6)
    p.add_argument('--radii',type=float,nargs='+',default=[2e-5,2e-4,2e-3])
    p.add_argument('--width-fractions',type=float,nargs='+',default=[0.,.01,.1])
    p.add_argument('--radial-orders',type=int,nargs='+',default=[24,32])
    p.add_argument('--core-orders',type=int,nargs='+',default=[12,20])
    p.add_argument('--ramp-fraction',type=float,default=0.)
    p.add_argument('--subcells',type=int,default=8)
    args = p.parse_args()
    if (not 0.<=args.ramp_fraction<.5 or args.subcells<2 or not args.q>0.
        or any(not 0.<r<.01 for r in args.radii) or any(not 0.<=g<=1. for g in args.width_fractions)
        or min(args.radial_orders+args.core_orders)<2):
        p.error('Invalid q, window, order or ramp controls.')
    fixture()
    print(json.dumps(run(**vars(args)),indent=2,allow_nan=False))
