"""Two equivalent delta-surface reductions of a frozen isotropic unit-residue test."""
import gzip
import json
from pathlib import Path
import functools

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import Q, intervals
from finite_mass_qed_archived_aperture_reference import difference, solve_window
from finite_mass_qed_shell_edge_slope_reference import SmoothFrozen


def integrate(P,B,m,Q,g,logs):
    g2=g*g
    geom=solve_window(P,B,m,Q,g2,logs)
    p0,p1=geom['endpoints']
    k0,k1=p0+Q,p1-Q
    Es=np.hypot(m,Q)
    if not k0<k1:
        raise ValueError('Kinematic chart does not invert the window.')
    def err(p,k):
        return difference(p,k,P.m,B.m,Es)+g2*(P.factor(p,logs)-B.factor(k,logs))
    def root(fun,l,r):
        return brentq(fun,l,r,rtol=1e-15,xtol=9e-15)
    def k_at(p): return root(lambda k:err(p,k),p-Q,p+Q)
    def p_at(k): return root(lambda p:err(p,k),k-Q,k+Q)

    @functools.lru_cache(maxsize=4096)
    def speed(p,parent):
        F=P if parent else B
        if logs and g:
            return F.pole.slope(p,g)
        E=np.hypot(p,F.m)
        return (p/E)*(1-g2*(F.photon+F.moment)/(2*np.pi**2*E**2))

    def internal_edges(F,L,R):
        q=np.r_[F.edges,F.pole.l,F.pole.r[-1]]
        return np.unique(q[(L<q)&(q<R)])

    ps=internal_edges(P,p0,p1)
    ks=internal_edges(B,k0,k1)
    pa=np.sort(np.r_[p0,p1,ps,[p_at(k) for k in ks]])
    ka=np.sort(np.r_[k0,k1,ks,[k_at(p) for p in ps]])
    if len(pa)!=len(ka):
        raise ValueError('The two chart representations use different support cells.')
    def run_in_chart(bounds,forward,n):
        t,w=roots_legendre(n)
        s=np.zeros((2,2))
        vmin=2.
        for l,r in zip(bounds[:-1],bounds[1:],strict=True):
            for a,v in zip((l+r)/2+(r-l)*t/2,(r-l)*w/2,strict=True):
                p,k=(a,k_at(a)) if forward else (p_at(a),a)
                va,vb=speed(p,True),speed(k,False)
                vmin=min(vmin,va,vb)
                if min(va,vb)<=.99:
                    raise ValueError('Group slope too small for either chart.')
                ep,eb=np.hypot(p,P.m),np.hypot(k,B.m)
                fac=p*k/((vb if forward else va)*ep*eb)
                # An optional denominator convention sensitivity, not a
                # prediction for the actual propagator residues.
                alt=ep*eb/(P.energy(p,g2,logs)*B.energy(k,g2,logs))
                f,b=P.f_at(p,True),B.f_at(k,True)
                s+=v*fac*np.outer([1.,alt],[f*(1-b),f-b])
        return s,vmin
    a,va=run_in_chart(pa,True,8)
    b,vb=run_in_chart(ka,False,12)
    defect=max(abs(a-b).ravel())
    if defect>2e-14:
        raise ValueError(f'Delta-surface chart mismatch {defect}.')
    if not g:
        # Independent Lorentz-invariant limit via an E_parent integrand.
        # In this limit the delta Jacobian k/(E_b v_b) is exactly 1.
        _,lo,hi,*_=intervals(P.m,m,B.m,Q)
        t,w=roots_legendre(14)
        c=np.hypot(pa,P.m)
        ref=np.zeros(2)
        for l,r in zip(c[:-1],c[1:],strict=True):
            ep=(l+r)/2+(r-l)*t/2
            pp=np.sqrt((ep-P.m)*(ep+P.m))
            k=np.sqrt((ep-Es-B.m)*(ep-Es+B.m))
            f=np.array([P.f_at(p,True) for p in pp])
            b=np.array([B.f_at(p,True) for p in k])
            ref+=(r-l)/2*np.array([w@(f*(1-b)),w@(f-b)])
        vac_err=max(abs(ref-a[0]))
        if abs(c[0]-lo)>2e-14 or abs(c[-1]-hi)>2e-14 or vac_err>2e-14:
            raise ValueError('The two charts fail the independent vacuum parent-energy integral.')
    else:
        vac_err=None
    return dict(bare_leg_integrals=a[0].tolist(),dressed_leg_integrals=a[1].tolist(),
                bare_vs_dressed_stimulus=a[1,1]/a[0,1]-1,different_chart_abs_defect=float(defect),
                vacuum_parent_energy_defect=vac_err,minimum_sampled_group_slope=min(va,vb),
                support_intervals=len(pa)-1,flat_stimulus=geom['limited_stimulus'])


def run():
    paths={g:Path(__file__).parent/f'results/early-pole-spectrum-{g}-m06-2026-11-26.json.gz'
           for g in ('fine','third')}
    items=[]
    for grid,p in paths.items():
        with gzip.open(p) as f:d=json.load(f)
        M=d['x_end']; b=M*d['mass_ratio']; m=d['checkpoints'][-1]['scalar_pole_over_T']
        for cut in (.25,.5):
            F=SmoothFrozen(d,'parent',M,cut,.05)
            B=SmoothFrozen(d,'fermion',b,cut,.05)
            vac=integrate(F,B,m,Q,0.,False)
            for g in (.005,.01):
                no=integrate(F,B,m,Q,g,False)
                yes=integrate(F,B,m,Q,g,True)
                diff=yes['bare_leg_integrals'][1]/no['bare_leg_integrals'][1]-1
                d2=yes['flat_stimulus']/no['flat_stimulus']-1
                items.append(dict(grid=grid,cut=cut,g=g,vacuum=vac,factorized=no,with_log=yes,
                                  log_weighted_stimulus_change=diff,log_flat_stimulus_change=d2,
                                  difference_of_relative_responses=diff-d2))
    return dict(scope='Frozen energy and geometry with explicit leg conventions; no matched collision amplitude.',
                rows=items)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
