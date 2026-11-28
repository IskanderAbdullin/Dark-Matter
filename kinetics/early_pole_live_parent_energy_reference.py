"""Continuity ledger with a live parent and configurable parent redistribution."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import X
from early_pole_parent_elastic_relaxation_reference import ElasticBath
from early_pole_parent_kinetic_relaxation_reference import KineticBath
from early_pole_resolved_parent_reference import ResolvedParent


def shell_mean_energy(d,x):
    qE=np.hypot(x,d.pq)
    return np.sum(d.pqw*qE,axis=1)/d.pw


def state_moments(d,x,y):
    a,b,P=np.split(y,(d.ns,d.iP))
    m0=1e-9*x
    m=d.d.gap(x,a)
    E=np.hypot(m,d.d.p)
    Ef=(d.daughter_energies(x) if hasattr(d,'daughter_energies')
        else .75*np.diff(d.d.fed**4)/np.diff(d.d.fed**3))
    Ep=shell_mean_energy(d,x)
    n=d.w*y
    U=(m*m-m0*m0)**2/(2*d.d.lam)
    daughter=n[:d.ns]@E-U+n[d.ns:d.iP]@Ef
    parent=n[d.iP:]@Ep
    return (np.r_[E,Ef,Ep],daughter,parent,m)


def rate(d,x,y):
    En,daughter,parent,m=state_moments(d,x,y)
    i,j,k,geom=d.overlaps(x,m)
    a,b,P=y[i],y[j],y[k]
    forward=geom*P*(1+a)*(1-b)
    reverse=geom*(1-P)*a*b
    R=forward-reverse
    # Zero physical energy imbalance in any exact continuum reaction;
    # this term isolates the first-moment residual of centroid updates.
    defect=(En[i]+En[j]-En[k])@R
    gross=En[k]@forward
    removed=En[k]@R
    m0=1e-9*x
    trace_scalar=m0*m0*((d.w[:d.ns]*y[:d.ns])@(1/En[:d.ns]))
    trace_p=x*x*np.sum(d.pqw/np.hypot(x,d.pq),axis=1)@y[d.iP:]
    heat=0.
    if isinstance(d,KineticBath):
        q,_,_=d.thermal_target(x,y[d.iP:])
        heat=d.rate*((d.pw*En[d.iP:])@(q-y[d.iP:]))
    elif isinstance(d,ElasticBath):
        q=d.elastic_target(x,y[d.iP:])
        heat=d.rate*((d.pw*En[d.iP:])@(q-y[d.iP:]))
    return np.array([defect,gross,removed,trace_p,trace_scalar,heat])


def ledger(ns,nf,nparent,mode,xend,mix):
    if mode=='none':
        d=ResolvedParent(ns,nf,nparent,0.,False)
    else:
        d=(KineticBath if mode=='fixed-T' else ElasticBath)(ns,nf,nparent,mix)
    init=np.r_[np.zeros(d.iP),d.parent_fd(X)]
    atol=np.r_[np.full(d.ns,1e-7),np.full(d.nf+d.np,8e-12)]
    sol=solve_ivp(d.rhs,(0.,np.log(xend/X)),init,method='BDF',jac=d.jac,
                  atol=atol,rtol=2e-8,max_step=.020,dense_output=True)
    assert sol.success,sol.message
    _,i_d,i_p,_=state_moments(d,X,init)
    _,f_d,f_p,m=state_moments(d,xend,sol.y[:,-1])
    output=[]
    for order in (2,4):
        q,w=roots_legendre(order)
        dt=np.diff(sol.t)/2
        eta=(sol.t[:-1,None]+dt[:,None]*(q+1)).ravel()
        we=(dt[:,None]*w).ravel()
        out=np.array([rate(d,X*np.exp(z),sol.sol(z)) for z in eta])
        output.append(we@out)
    F=output[-1]
    discrepancy=np.max(abs(F-output[0]))/f_d
    assert discrepancy<2e-5, (F,output[0])
    dp=f_p-i_p
    dd=f_d-i_d
    defect,gross,removed,tp,ts,heat=F
    closure_all=(dp+dd-(defect+tp+ts+heat))/dd
    closure_daughters=(dd-(removed+defect+ts))/dd
    closure_parent=(dp-(-removed+tp+heat))/dd
    for e in (closure_all,closure_daughters,closure_parent):
        assert abs(e)<5e-5,(e,output)
    return dict(model=mode,rate=mix if mode!='none' else 0,
                cells=[d.ns,d.nf,d.np],xend=xend,
                energy_unit='rho/Tref^4; Tref proportional to 1/a',
                initial_parent_energy=i_p,final_parent_energy=f_p,
                final_daughter_energy=f_d,final_scalar_pole_over_T=m,
                signed_decay_centroid_defect=float(defect),
                gross_parent_decay_moment=float(gross),net_parent_decay_moment=float(removed),
                integrated_parent_trace=float(tp),integrated_scalar_trace=float(ts),
                integrated_heat_from_fixed_temperature_bath=float(heat),
                centroid_defect_fraction_of_final_daughter_energy=float(defect/dd),
                signed_heat_fraction_of_final_daughter_energy=float(heat/dd),
                absolute_continuity_errors_over_final_daughter_energy=dict(
                    total=float(abs(closure_all)),daughters=float(abs(closure_daughters)),
                    parent=float(abs(closure_parent))),
                two_order_quadrature_disagreement_over_daughter_energy=float(discrepancy))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--mode',choices=['none','fixed-T','isolated'],default='none')
    p.add_argument('--rate',type=float,default=100.)
    p.add_argument('--xend',type=float,default=.05)
    a=p.parse_args()
    assert a.rate>=0 and X<a.xend<=.05
    print(json.dumps(ledger(a.scalar,a.fermion,a.parent,a.mode,a.xend,a.rate),
                     indent=2,allow_nan=False))
