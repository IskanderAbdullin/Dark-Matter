"""Evaluate momentum-local withdrawals demanded of the held FD parent bath."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import G, KAPPA, X
from early_pole_paired_hartree_matrix_reference import DirectDriver


def stable_integral(ll, dr):
    """Integral of the unit-temperature FD distribution over a short interval."""
    return np.log1p(np.exp(-ll)*(-np.expm1(-dr))/(1+np.exp(-ll-dr)))


def demands(d, x, u, b, parent_edges):
    m = d.gap(x, u)
    E = np.hypot(m, d.p)
    Delta = (x*x-m*m)/2
    lo = (x*x*E*E+Delta*Delta)/(((x*x+m*m)/2)*E+d.p*Delta)
    hi = lo+2*d.p*Delta/(m*m)
    assert x > m
    parent = np.hypot(x, parent_edges)
    daughter = E[:, None] + d.ef
    # The sorted union changes both piecewise daughter f and parent cell only
    # at edges. Clip each elementary subinterval to the scalar's aperture.
    cut = np.sort(np.concatenate((np.broadcast_to(parent, (d.n,len(parent))),
                                  daughter), axis=1),axis=1)
    lower = np.maximum(cut[:, :-1],lo[:,None])
    upper = np.minimum(cut[:, 1:],hi[:,None])
    dr = np.maximum(upper-lower,0.)
    mid = lower + dr/2
    q = mid - E[:, None]
    j = np.searchsorted(d.ef, q, side='right')-1
    k = np.searchsorted(parent, mid, side='right')-1
    valid = (j >= 0) & (j < len(b)) & (k >= 0) & (k < len(parent)-1)
    fp = stable_integral(lower, dr)
    b_at = b[np.clip(j,0,len(b)-1)]
    row = d.w * KAPPA * G * x * (x*x-m*m)/(E*d.p)
    forward = (row*(1+u))[:,None]*fp*(1-b_at)
    reverse = (row*u)[:,None]*(dr-fp)*b_at
    out = np.zeros((2,len(parent)-1))
    np.add.at(out[0],k[valid],forward[valid])
    np.add.at(out[1],k[valid],reverse[valid])
    return out


def integrate(ns,nf,xend,bins):
    d = DirectDriver(ns,nf,7.65175e-8,0.,1e-9)
    eta = np.log(xend/X)
    atol = np.r_[np.full(d.n,5e-8),np.full(len(d.wf),6e-12)]
    sol = solve_ivp(d.rhs,(0.,eta),np.zeros(len(atol)),method='BDF',jac=d.jac,
                    atol=atol,rtol=1e-8,max_step=.02,dense_output=True)
    assert sol.success,sol.message
    edges = np.r_[0.,np.geomspace(1e-3,8.,bins),32.,60.]
    # Gauss quadrature separately across the solver's adaptive segments.
    # Save cumulative demand at the endpoint and at an earlier comparable x.
    all_checks = []
    for order in (2,4):
        zz,ww = roots_legendre(order)
        h = np.diff(sol.t)/2
        t = sol.t[:-1,None]+h[:,None]*(1+zz)
        W = h[:,None]*ww
        fold = np.zeros((len(sol.t)-1,2,len(edges)-1))
        for i,(ts,ws) in enumerate(zip(t,W,strict=True)):
            for z,w in zip(ts,ws,strict=True):
                v=sol.sol(z)
                fold[i] += w*demands(d,X*np.exp(z),v[:d.n],v[d.n:],edges)
        F = fold.cumsum(axis=0)
        all_checks.append(F)
    F = all_checks[-1]
    net = F[:,0]-F[:,1]
    assert np.max(np.abs(F[-1]-all_checks[0][-1])) < 2e-6*np.max(F[-1])
    Q,W = roots_legendre(64)
    q = edges[:-1,None]+np.diff(edges)[:,None]*(Q+1)/2
    w = 4*np.diff(edges)[:,None]*W*q*q/(4*np.pi**2)
    def fd_inventory(x):
        return np.sum(w/(1+np.exp(np.hypot(q,x))),axis=1)
    n0 = fd_inventory(X)
    nend = fd_inventory(xend)
    assert abs(net[-1].sum()/(d.w@sol.y[:d.n,-1])-1)<3e-5
    mask=n0>1e-14
    k=int(np.argmax(np.where(mask,net[-1]/n0,-np.inf)))
    peak=int(np.argmax(np.max(np.where(mask,net/n0,-np.inf),axis=1)))
    return dict(method='held FD parent; signed local demand, not actual depleted spectrum',
                parent_comoving_momentum_edges_over_T=edges.tolist(),
                initial_parent_n_over_T3=n0.tolist(),
                equilibrium_endpoint_parent_n_over_T3=nend.tolist(),
                integrated_forward_parent_n_over_T3=F[-1,0].tolist(),
                integrated_inverse_parent_n_over_T3=F[-1,1].tolist(),
                global_requested_fraction=float(net[-1].sum()/n0.sum()),
                strongest_endpoint_cell=[float(edges[k]),float(edges[k+1])],
                strongest_endpoint_cell_requested_fraction=float(net[-1,k]/n0[k]),
                strongest_cell_time_peak_requested_fraction=float(np.max(np.where(mask,net/n0,-np.inf))),
                peak_x=float(X*np.exp(sol.t[1:][peak])),
                n_cells_exceeding_inventory_at_endpoint=int(np.count_nonzero(mask & (net[-1]>n0))),
                inventory_above_noise_floor_over_T3=1e-14,
                worst_energy_support_omission_test_relative=float(abs(net[-1].sum()/(d.w@sol.y[:d.n,-1])-1)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=65)
    p.add_argument('--fermion',type=int,default=90)
    p.add_argument('--bins',type=int,default=36)
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    print(json.dumps(integrate(a.scalar,a.fermion,a.xend,a.bins),indent=2,allow_nan=False))
