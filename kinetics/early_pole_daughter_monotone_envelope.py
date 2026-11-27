"""Monotone per-cell daughter occupations and a cutoff-mixture dual."""
import argparse
import json
from math import fsum

import numpy as np
from scipy.optimize import linprog
from scipy.special import expit,roots_legendre

from early_pole_daughter_frozen_moment_envelope import cell_objects
from early_pole_daughter_local_fraction_curvature_attribution import initial


def upper_lower(m,obj,physical_measure,scale,sample,relax=0.):
    A=np.c_[np.zeros(len(m)),np.cumsum(physical_measure,axis=1)]/scale[0]
    R=np.r_[0.,np.cumsum(obj)]/scale[1]
    eq=np.vstack((np.ones(len(R)),A))
    target=np.r_[1.,m/scale[0]]
    vals=[]
    for s in (1.,-1.):
        opt=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)
        if relax:
            B=np.r_[target[1:]+relax,-target[1:]+relax]
            lp=linprog(s*R,A_eq=eq[:1],b_eq=target[:1],A_ub=np.vstack((A,-A)),b_ub=B,
                       bounds=(0.,1.),method='highs',options=opt)
        else:
            lp=linprog(s*R,A_eq=eq,b_eq=target,bounds=(0.,1.),method='highs',options=opt)
        if not lp.success:
            raise ValueError(lp.message)
        if relax:
            mu=np.minimum(lp.ineqlin.marginals,0.)
            k=len(m)
            lam=(mu[:k]-mu[k:])*scale[1]/scale[0]
            constant=float(mu@B*scale[1])
        else:
            lam=lp.eqlin.marginals[1:]*scale[1]/scale[0]
            constant=float(lam@m)
        limits=[]
        for order in (32,64):
            w,AA,k,_=sample(order)
            t=s*k-np.einsum('i,iab->ab',lam,AA)
            # Cutoff cost H(t)-lambda*M(t) is an integrated density.
            # Within a bin, lower-bound it using its full negative/positive mass
            # from each endpoint. Take the best endpoint bound in each bin.
            k0=np.r_[0.,np.cumsum(np.sum(w*t,axis=1))]
            lower=np.maximum(k0[:-1]+np.sum(w*np.minimum(t,0.),axis=1),
                             k0[1:]-np.sum(w*np.maximum(t,0.),axis=1))
            limits.append(s*float(constant+min(lower)))
        vals.append(dict(finite_lp_objective=float(s*lp.fun*scale[1]),
                         numerical_continuous_bound=limits[-1],
                         dual_order_32_to_64_delta=abs(limits[-1]-limits[0]),
                         finite_normalized_moment_residual=float(max(abs(eq@lp.x-target)))))
    return vals


def run(archive,sub,gauge):
    C,n,(y0,_,y1),_,depth=initial(archive)
    M=C.M
    e,f=M.weights(y0,True)[:2]
    pole=C.scalar(n)
    aa=np.asarray(M.edge_sets[2])
    cells=M.r[2].start+np.flatnonzero((aa[:-1]>=1.4138888888888888-2e-14)&(aa[1:]<=1.4666666666666666+2e-14))
    if len(cells)!=2*depth:
        raise ValueError('Wrong physical daughter interval.')
    data=[]
    u,g=roots_legendre(64)
    for j in cells:
        B=M.cells[j]
        sample,A,W,den,K,_=cell_objects(C,int(j),e,f,sub,pole,n,gauge,min_half_width=0.)
        w,features,k,p=sample(32)
        scale=(B.hi-B.lo)/2
        P=(B.lo+B.hi)/2+scale*u
        full=g*scale*2*P*P/np.pi**2
        q=float(y0[j,0]/full.sum())
        err=abs(y0[j,1]-q*(full@B.basis(P)))
        if err>8e-17 or abs(B.fit(y0[j])[0][1])>1e-12:
            raise ValueError('Native first two moments are not the constant-occupation moments.')
        z=B.fit(y1[j])[0]
        if z[1]>=0.:
            raise ValueError('The local reconstructed occupation is not decreasing in this fixture.')
        fp=expit(z[0]+z[1]*B.basis(p))
        target=np.r_[y1[j],np.sum(w*features[2]*fp)] if gauge else y1[j]
        lp=upper_lower(target,W,A,(den,K),sample)
        m0=np.r_[y0[j],q*np.sum(w*features[2])] if gauge else y0[j]
        rob=upper_lower(m0,W,A,(den,K),sample,relax=1e-6)
        native=float(q*np.sum(w*k))
        fit=float(np.sum(w*k*fp))
        a,b=lp
        if not (a['numerical_continuous_bound']-1e-9<=fit<=b['numerical_continuous_bound']+1e-9):
            raise ValueError('A smooth monotone fitted control is outside its proposed dual envelope.')
        data.append(dict(cell=int(j),flat_native_energy_moment_defect=float(err),
                         local_monotone_fitted_influence=fit,native_flat_influence=native,
                         local_cutoff_min=a,local_cutoff_max=b,
                         native_relaxed_min=rob[0],native_relaxed_max=rob[1]))
    lower=fsum(v['native_flat_influence']-v['local_cutoff_max']['numerical_continuous_bound'] for v in data)
    upper=fsum(v['native_flat_influence']-v['local_cutoff_min']['numerical_continuous_bound'] for v in data)
    return dict(scope='Monotone nonincreasing within each fixed daughter cell; jumps between cells permitted; frozen packet geometry.',
                hard_depth=depth,cold_modes=len(n),subdivision=sub,fixed_gauge_moment=gauge,
                independent_fitted_pair=fsum(v['native_flat_influence']-v['local_monotone_fitted_influence'] for v in data),
                numerical_outer_paired_interval=[lower,upper],native_normalized_moment_relaxation=1e-6,
                relaxed_native_numerical_outer_interval=[
                    fsum(v['native_relaxed_min']['numerical_continuous_bound']-v['local_cutoff_max']['numerical_continuous_bound'] for v in data),
                    fsum(v['native_relaxed_max']['numerical_continuous_bound']-v['local_cutoff_min']['numerical_continuous_bound'] for v in data)],
                cells=data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archives',nargs='+')
    p.add_argument('--subdivision',type=int,default=12)
    p.add_argument('--fix-gauge',action='store_true')
    a=p.parse_args()
    if a.subdivision<1:
        p.error('Require positive cutoff subdivisions.')
    print(json.dumps([run(f,a.subdivision,a.fix_gauge) for f in a.archives],
                     indent=2,allow_nan=False))
