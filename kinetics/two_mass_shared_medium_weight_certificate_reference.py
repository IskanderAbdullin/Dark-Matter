"""Fixed-extent weight-box bounds and contraction for a quadratic medium."""
import itertools
import json
from pathlib import Path

import numpy as np

from two_mass_kinetic_offset_deposition_reference import DT, event


def certificate(old_gradient, hessian, a, d, lo, hi, xi, capacity, n, fermi):
    """Directions v_j(w_j)=a_j+w_j*d_j; v_j and d_j are stored in rows."""
    left, right = a + lo[:, None]*d, a + hi[:, None]*d
    ends = np.stack((left, right))
    # These are extrema over independent weight coordinates at fixed extents.
    # Only the k-th direction depends on the k-th coordinate.
    def bounds(probes):
        coeff = (xi/2)*(probes @ hessian @ ends.swapaxes(1, 2))
        return (probes@old_gradient + coeff.min(axis=0).sum(axis=-1),
                probes@old_gradient + coeff.max(axis=0).sum(axis=-1))
    lmin, lmax = bounds(left)
    umin, umax = bounds(right)
    minimum_slope, maximum_slope = bounds(d)
    cross = np.asarray([(xi/2)*(probes @ hessian @ ends.swapaxes(1,2)) for probes in ends])
    radius = np.max(np.abs(cross),axis=1).sum(axis=-1)
    margin = np.minimum(-left@old_gradient,right@old_gradient)
    symmetric = min(margin-radius.max(axis=0))
    map_margin = min((-lmax).min(),umin.min())
    # If one event cannot reach zero anywhere, there cannot be a common root.
    exclusion = max(lmin.max(), (-umax).max()) if min(minimum_slope)>0 else DT(0)
    # Bound implicit derivative of the frozen-gradient root map.
    # Evaluate the output direction at endpoints to bound its coefficient.
    numer = np.max(np.abs(ends @ hessian @ d.T),axis=0)
    contractions = np.full_like(minimum_slope,np.inf)
    np.divide(np.sum(numer*np.abs(xi/2),axis=1),minimum_slope,
              where=minimum_slope>0,out=contractions)
    increment = xi[None,:,None]*ends
    state_lo = n + increment.min(axis=0).sum(axis=0)
    state_hi = n + increment.max(axis=0).sum(axis=0)
    occupancy_margin = min(state_lo.min(), (capacity-state_hi)[fermi].min())
    return dict(lmin=lmin,lmax=lmax,umin=umin,umax=umax,min_slope=minimum_slope,
                max_slope=maximum_slope,symmetric_margin=symmetric,
                map_margin=map_margin,exclusion_margin=exclusion,
                lipschitz=contractions,occupancy_margin=occupancy_margin,
                state_lo=state_lo,state_hi=state_hi)


def run():
    e,m,b,P,g=map(DT,(.046,.002,.007,.36,.012))
    ep=np.hypot(e+m+b,P)
    es,fd=np.array([event(e,m,b,P,x)[1]+np.array([m,b],dtype=DT)
                    for x in (-.38,-.31)]).T
    bare=np.array([ep,min(fd)*DT(.94),max(fd)*DT(1.06),*es],dtype=DT)
    gamma=g*g/DT(np.pi)**4
    B=np.zeros((2,5),dtype=DT)
    B[0,0]=DT(np.pi)**2/(2*ep)
    B[1,1:3]=DT(np.pi)**2/(2*bare[1:3])
    J0=np.asarray([.62,.32],dtype=DT)+DT(.27)
    H=gamma*B.T@B
    cap=np.array([.9,.8,.9,1.,1.],dtype=DT)
    n=cap*np.array([.16,.29,.31,.005,1.1],dtype=DT)
    fermi=np.asarray([True,True,True,False,False])
    a=np.array([[-1.,1.,0.,1.,0.],[-1.,1.,0.,0.,1.]],dtype=DT)
    d=np.broadcast_to(np.array([0.,-1.,1.,0.,0.],dtype=DT),a.shape)

    def grad(z):
        return bare+gamma*B.T@(J0+B@z)

    gn=grad(n)
    center=-(a@gn)/(d@gn)
    data=Path(__file__).parent/'results/two-mass-shared-medium-event-weights-2026-11-27.json'
    result=[]
    for row in json.loads(data.read_text())['rows']:
        xi=np.asarray(row['extents'],dtype=DT)
        for width in (4e-5,6e-5,None):
            lo,hi=(np.zeros(2,dtype=DT),np.ones(2,dtype=DT)) if width is None else (center-width,center+width)
            cert=certificate(gn,H,a,d,lo,hi,xi,cap,n,fermi)
            left,right=a+lo[:,None]*d,a+hi[:,None]*d
            def state(w):
                return n+xi@(a+w[:,None]*d)
            # Independent evaluation at all corners using the unsplit gradient.
            samples=np.asarray([state(np.asarray(w)) for w in itertools.product(*zip(lo,hi,strict=True))])
            grids=np.asarray([grad((n+y)/2) for y in samples])
            comparisons=[(cert['state_lo'],samples.min(axis=0)),(cert['state_hi'],samples.max(axis=0))]
            for probes,k0,k1 in ((left,'lmin','lmax'),(right,'umin','umax'),(d,'min_slope','max_slope')):
                val=grids @ probes.T
                comparisons.extend(((cert[k0],val.min(axis=0)),(cert[k1],val.max(axis=0))))
            mismatch=max(np.max(np.abs(x-y)) for x,y in comparisons)
            if mismatch>1e-16 or cert['occupancy_margin']<=0 or min(cert['min_slope'])<=0:
                raise ValueError('Corner-extrema/occupation/slope test failed.')
            # Require ample guard relative to tested floating-point error.
            q=max(cert['lipschitz'])
            certified=bool(cert['map_margin']>1e-14 and q<1)
            excluded=bool(cert['exclusion_margin']>1e-14)
            iterations=None
            agrees=None
            if certified:
                ww=(lo+hi)/2
                for i in range(1,31):
                    gm=grad((n+state(ww))/2)
                    nw=-(a@gm)/(d@gm)
                    if np.any(nw<lo) or np.any(nw>hi) or excluded:
                        raise ValueError('Certified weight map left its box.')
                    move=max(abs(nw-ww))
                    ww=nw
                    if move<1e-17:
                        iterations=i
                        break
                if iterations is None:
                    raise ValueError('Contraction failed to converge.')
                agrees=float(max(abs(ww-np.asarray(row['weights'],dtype=DT))))
                if agrees>3e-15 or max(abs((a+ww[:,None]*d)@grad((n+state(ww))/2)))>1e-16:
                    raise ValueError('Picard disagrees with coupled benchmark.')
            result.append(dict(step=row['step'],half_width=width,
                               fixed_extent_common_root_certified=certified,
                               fixed_extent_common_root_excluded=excluded,
                               signed_map_margin=float(cert['map_margin']),
                               exclusion_margin=float(cert['exclusion_margin']),
                               symmetric_margin=float(cert['symmetric_margin']),
                               maximum_contraction_bound=float(q),
                               minimum_occupation_margin=float(cert['occupancy_margin']),
                               corner_extrema_max_abs_error=float(mismatch),
                               picard_iterations=iterations,
                               picard_weight_error_vs_coupled_reference=agrees))
    if not any(r['fixed_extent_common_root_certified'] and r['symmetric_margin']<0 for r in result):
        raise ValueError('Expected a cancellation-resolved control.')
    if not any(r['fixed_extent_common_root_excluded'] for r in result):
        raise ValueError('Expected an excluded narrow box.')
    return dict(scope='Two fixed-extent events in a prescribed quadratic medium; no thermal matching.',
                rows=result)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
