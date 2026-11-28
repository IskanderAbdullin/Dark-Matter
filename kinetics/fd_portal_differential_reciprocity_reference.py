"""Differential FD number/heat response away from common-bath detailed balance."""
import argparse
import json
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit

from positive_fd_portal_quadrature_reference import PositivePortal, rule, leptons


def integrate(p,Ts,alpha,beta,*,invariant_window=None):
    """At fixed Ts, coordinates are 2 eta and 1 - Ts / Td."""
    if not Ts>0. or not beta<1. or not np.all(np.isfinite([Ts,alpha,beta])):
        raise ValueError('Need finite positive local temperature and admissible force')
    Td=Ts/(1-beta)
    E,q,p3,W,vi2=p.geometry(max(Ts,Td))
    if invariant_window is not None:
        a,b=invariant_window
        if not 0.<=a<b:
            raise ValueError('Need an ordered invariant-energy window')
        W=np.where((E[:,None]>=a)&(E[:,None]<b),W,0.)
    u,w=rule(p.orders[2])
    u=(u+1)/2
    w=w/2
    w2=w*(3*u*u-1)/2
    pair=q/Ts
    dq=p3[:,:,None]*np.sqrt(vi2)[:,None,None]*u/(2*Ts)
    e1=pair[:,:,None]/2+dq
    e2=pair[:,:,None]/2-dq
    f1=expit(alpha/2-(1-beta)*e1)
    f2=expit(alpha/2-(1-beta)*e2)
    D=(1-f1)*(1-f2)
    da=-.5*(f1+f2)*D
    db=-(e1*f1+e2*f2)*D
    anti=.5*(e1-e2)*(f2-f1)*D
    if np.min(anti)<0. or np.max(abs(db-pair[:,:,None]*da-anti))>8e-14:
        raise ArithmeticError('Exact sign of the angular blocking derivative failed')
    dark=[(a@w,a@w2) for a in (D,da,db,anti)]
    light=np.zeros((4,)+q.shape)
    s=E*E
    for m in leptons:
        ok=E>2*m
        vf=1-4*m*m/s[ok]
        v=vi2[ok]
        d=p3[ok,:,None]*np.sqrt(vf)[:,None,None]*u/(2*Ts)
        g=expit(pair[ok,:,None]/2+d)*expit(pair[ok,:,None]/2-d)
        pref=((1+2*m*m/s[ok])*np.sqrt(vf))[:,None]
        f=((2/3)*v*vf/(3-v-vf+v*vf/3))[:,None]
        if np.any(f<0.) or np.any(f>1.):
            raise ArithmeticError('Nonnegative spin-correlated angular kernel required')
        c0,c2=g@w,g@w2
        for j,(v0,v2) in enumerate(dark):
            light[j,ok]+=pref*(v0[ok]*c0+f*v2[ok]*c2)
    I,J,K,L=light*W
    d=alpha+beta*pair
    A=np.exp(alpha-(1-beta)*pair)
    B=np.exp(-pair)
    diff=np.exp(np.maximum(alpha-(1-beta)*pair,-pair))*np.sign(d)*(-np.expm1(-abs(d)))
    # At a fixed physical node only the population changes with
    # alpha and beta; the quadrature shape is not differentiated.
    C=I*diff
    f=np.array([np.sum(C),np.sum(C*pair),np.sum(C*d)])
    ac=I*A+J*diff
    bc=I*A*pair+K*diff
    D=np.array([[np.sum(ac),np.sum(bc)],[np.sum(pair*ac),np.sum(pair*bc)]])
    antisymmetric=np.sum(L*diff)
    if np.min(D)<=0. or f[2]<0. or abs(D[0,1]-D[1,0]-antisymmetric)>3e-14*np.max(D):
        raise ArithmeticError('Derivative positivity, entropy, or direct antisymmetry failed')
    return f,D,antisymmetric,np.array([np.sum(L*B),np.sum(L*B*pair)])


def zeros(p,T,b):
    def node(alpha):
        f,D,A,_=integrate(p,T,alpha,b)
        return f,D,A
    def root(target):
        def fn(x):
            f,D,A=node(x)
            return (f[0],f[1],A)[target]
        u=brentq(fn,-15.,15.,xtol=2e-13)
        f,D,A=node(u)
        den=np.sqrt(D[0,1]*D[1,0])
        if abs((f[0],f[1],A)[target])>2e-12*den:
            raise ArithmeticError('Selected differential or flow zero did not close')
        return u,A/den
    return [root(i) for i in range(3)]


def run():
    items=[(.05,0.,0.),(.3,0.,0.),(.6882487043,0.,0.),
           (.05,.6,.16),(.05,-.6,-.16),(.3,-1.5,.16),(.3,1.5,-.16),
           (.6882487043,-1.8,.16),(.6882487043,1.8,-.16)]
    a,b=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for T,u,v in items:
        lo,A,aa,_=integrate(a,T,u,v)
        f,D,anti,_=integrate(b,T,u,v)
        e=T/(1-v)
        q=b.moments(T,e,u/2)
        raw=q[:3]/[1.,T,1.]
        scale=np.max(abs(D))
        if np.max(abs(f-raw))>3e-13*max(np.max(abs(raw)),scale) or np.max(abs(A-D))/scale>2e-6:
            raise ArithmeticError('Analytic derivative did not close with the independent original FD moments or radial orders')
        def shifted(i,h):
            x=u+h*(i==0)
            y=v+h*(i==1)
            return b.moments(T,T/(1-y),x/2)[:3]/[1.,T,1.]
        step=8e-5
        N=np.array([(shifted(i,-2*step)-8*shifted(i,-step)+8*shifted(i,step)-shifted(i,2*step))/(12*step)
                    for i in (0,1)]).T
        derivative_entropy=np.array([f[0],f[1]])+np.array([u,v])@D
        error=np.max(abs(np.vstack((D,derivative_entropy))-N))/max(scale,np.max(abs(derivative_entropy)))
        if error>3e-7:
            raise ArithmeticError(('Independent finite differences did not match',T,u,v,float(error)))
        anti_scaled=anti/np.sqrt(D[0,1]*D[1,0])
        if u==v==0.:
            known=np.array([[q[3],q[4]/T],[q[4]/T,q[5]/(T*T)]])
            if anti!=0. or np.max(abs(D/known-1.))>2e-13 or np.linalg.det(D)<=0.:
                raise ArithmeticError('Common-bath Onsager coefficients did not match')
        elif u*v>=0. and np.sign(anti)!=np.sign(u+v):
            raise ArithmeticError('One-signed affinities contradict the strict direct antisymmetry')
        rows.append(dict(T_SM=T,T_D=e,coordinates=[u,v],number_heat_and_entropy=f.tolist(),
                         analytic_number_heat_jacobian=D.tolist(),antisymmetry_fraction=anti_scaled,
                         determinant_over_geometric_row_product=float(np.linalg.det(D)/(D[0,0]*D[1,1])),
                         symmetric_part_eigenvalues_normalized=(np.linalg.eigvalsh((D+D.T)/2)/scale).tolist(),
                         max_finite_difference_error=float(error),max_independent_order_disagreement=float(np.max(abs(D-A))/scale)))
    scans=[]
    for T in (.05,.30,.6882487043):
        _,D,_,W=integrate(b,T,0.,0.)
        slopes=np.array([-D[0,1]/D[0,0],-D[1,1]/D[1,0],-W[1]/W[0]])
        if np.min(W)<=0.:
            raise ArithmeticError('Need positive derivative weights')
        delta=3e-4
        eps=np.array([zeros(b,T,v) for v in (-delta,delta)])
        d=(eps[1,:,0]-eps[0,:,0])/(2*delta)
        if np.max(abs(d/slopes-1.))>3e-6:
            raise ArithmeticError('Independent nonlinear roots failed the common linear prediction')
        separated=[]
        for v in (-.12,.12):
            x,y=np.asarray(zeros(a,T,v)),np.asarray(zeros(b,T,v))
            if np.max(abs(x-y))>3e-6:
                raise ArithmeticError('Independent meshes do not agree on the selected zeros')
            rate,_,a0,_=integrate(b,T,y[2,0],v)
            if (not rate[2]>0. or np.min(abs(y[2,0]-y[:2,0]))<1e-5
                    or abs(rate[2]-(y[2,0]*rate[0]+v*rate[1]))>2e-13*np.max(abs(rate))):
                raise ArithmeticError('The symmetric differential state did not retain distinct positive entropy')
            separated.append(dict(beta=v,alpha_zeros_number_heat_reciprocity=y[:,0].tolist(),
                                  off_diagonal_asymmetry_fraction_at_zeros=y[:,1].tolist(),
                                  reciprocal_state_flows_and_entropy=rate.tolist()))
        scans.append(dict(T_SM=T,common_limit_alpha_derivative_in_beta=slopes.tolist(),
                          independent_root_derivative=d.tolist(),finite_displacements=separated))
    def common_slopes(q,T):
        _,D,_,W=integrate(q,T,0.,0.)
        return np.array([-D[0,1]/D[0,0],-D[1,1]/D[1,0],-W[1]/W[0]])
    Ts=np.linspace(.30,.6882487043,17)
    mats=[[common_slopes(q,T) for T in Ts] for q in (a,b)]
    transition=[]
    for j,label in enumerate(('number','heat')):
        keys=[]
        for q,w in zip((a,b),mats,strict=True):
            w=np.asarray(w)
            picks=np.flatnonzero(np.diff(np.sign(w[:,2]-w[:,j]))!=0)
            brackets=[]
            for i in picks:
                def difference(T):
                    z=common_slopes(q,T)
                    return z[2]-z[j]
                hit=brentq(difference,Ts[i],Ts[i+1],xtol=2e-13)
                brackets.append(dict(bracket=Ts[i:i+2].tolist(),visible_temperature=hit,
                                     slopes=common_slopes(q,hit).tolist()))
            if not brackets:
                raise ArithmeticError('No local root found inside the sampled interval')
            keys.append(brackets)
        if len(keys[0])!=len(keys[1]) or max(abs(x['visible_temperature']-y['visible_temperature'])
                                            for x,y in zip(*keys,strict=True))>3e-7:
            raise ArithmeticError('The common-limit local slope intersections did not resolve')
        transition.append(dict(shared_slope_with=label,independently_bracketed_entries=keys[1]))
    return dict(rows=rows,separate_differential_nulls=scans,common_limit_sampled_slope_intersections=transition,forces=['2 eta_D','1-T_SM/T_D'],flows=['F','Q/T_SM'],
                scope='Direct charged-lepton vacuum FD benchmark at fixed visible temperature. Away from detailed balance the derivative includes state-dependent blocking. This is not a thermodynamic trajectory, abundance, or decoupling diagnosis.')


def intermediate():
    T=.64135
    endpoints=[]
    q=(PositivePortal((13,36,26)),PositivePortal((17,44,32)))
    for b in (-.002,.002):
        row=[]
        for p in q:
            u=np.array([a[0] for a in zeros(p,T,b)])
            y,D,skew,_=integrate(p,T,u[2],b)
            if (not min(u[:2])<u[2]<max(u[:2]) or not y[0]*y[1]<0. or y[2]<=0.
                    or abs(skew)>2e-12*np.max(D)):
                raise ArithmeticError('Independent intermediate differential zero not resolved')
            row.append(dict(orders=p.orders,alpha_zeros=u.tolist(),
                            flows_number_heat_over_T_entropy=y.tolist()))
        A,B=[np.asarray(p['flows_number_heat_over_T_entropy']) for p in row]
        if np.max(abs(A/B-1.))>2e-6:
            raise ArithmeticError('Intermediate actual fluxes fail two orders')
        endpoints.append(dict(beta=b,by_independent_order=row))
    return dict(visible_temperature=T,both_signs=endpoints,
                scope='Direct small finite displacements between the two sampled common-slope intersections. The differential cross terms balance while the two finite fluxes have opposite signs and entropy is positive.')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--intermediate',action='store_true')
    print(json.dumps(intermediate() if a.parse_args().intermediate else run(),indent=2,allow_nan=False))
