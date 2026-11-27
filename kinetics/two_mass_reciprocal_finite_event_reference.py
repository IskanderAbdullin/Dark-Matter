"""Finite event balance for a specified reciprocal quadratic medium proxy."""
import json

import numpy as np
from scipy.optimize import brentq

from early_pole_reciprocal_gauge_energy_functional_reference import setup,state
from two_mass_kinetic_offset_deposition_reference import DT,event
from two_mass_entropy_extent_reference import entropy
from two_mass_discrete_entropy_gradient_reference import averaged_alphas


def model(name,e,m,b,momentum,angle,g):
    e,m,b,momentum,angle,g=map(DT,(e,m,b,momentum,angle,g))
    x=m+b+e
    _,K,avail=event(e,m,b,momentum,angle)
    ep=x+avail-e
    es=m+K[0]
    ed=b+K[1]
    if abs(ep-es-ed)>3e-16*ep or ed*DT(.96)<=b:
        raise ValueError("Vacuum event or lower receiving representative failed the kinematic check.")
    # Two receiving fermion representatives surrounding the exact bare
    # two-body target. The scalar and parent are held representatives.
    ebase=np.array([ep,ed*.96,ed*1.04,es],dtype=DT)
    q,w,f,d=setup()
    J_ref=state(1.,np.array([0.,float(x),float(b)]),float(g),q,w,f,d)[-1]
    A=g*g/DT(np.pi)**4
    B=np.zeros((2,4),dtype=DT)
    B[0,0]=DT(2)*DT(np.pi)**2/(DT(d[1])*ep)
    B[1,1:3]=DT(2)*DT(np.pi)**2/(DT(d[2])*ebase[1:3])
    n=np.array([.045,.025,.035,.018],dtype=DT)
    J0=np.asarray(J_ref[1:],dtype=DT)
    photon=DT(J_ref[0])
    js=lambda z:J0+B@z
    def total(z):
        J=js(z)
        return ebase@z+A*(photon*J.sum()+J@J/2)
    def gradient(z):
        return ebase+A*B.T@(photon+js(z))
    H=A*B.T@B
    def direction(w):
        return np.array([-1,1-DT(w),DT(w),1],dtype=DT)
    def residual(w,dx,origin):
        v=direction(w)
        return gradient(origin+dx*v/2)@v
    def weight(dx,origin=n):
        lo,hi=[residual(z,dx,origin) for z in (0.,1.)]
        if not lo<0<hi:
            raise ValueError("The finite event is outside the positive energy bracket.")
        return brentq(lambda z:float(residual(z,dx,origin)),0.,1.,xtol=5e-15,rtol=1e-15)
    base=weight(DT(0))
    trace_flags=[]
    for dx in map(DT,(-.012,-.003,.003,.012)):
        settled=weight(dx)
        v=direction(settled)
        n1=n+dx*v
        if min(n1)<=0 or max(n1)>=1 or np.max(abs(np.array([1,1,1,0])@v))>1e-18 or abs((np.array([1,0,0,1])@v))>1e-18:
            raise ValueError("The finite event lost a conserved count or the manufactured physical box.")
        eta0=direction(base)
        old=total(n+dx*eta0)-total(n)
        theory=dx*dx*(eta0@H@eta0)/2
        fixed=total(n1)-total(n)
        # Absolute and normalized finite energy checks, plus independent
        # verification that the missing drift is quadratic.
        if abs(fixed)>3e-16*total(n) or theory<=0 or abs(old/theory-1)>1.3e-8:
            raise ValueError("Frozen and exact medium ledgers did not separate as predicted.")
        if abs(settled-base)<=1e-13:
            raise ValueError("Medium finite-step weight correction was not resolved.")
        # Reverse at the exact new state. The same midpoint is reached
        # when swapping endpoints, so solve it independently there.
        inverse=lambda s:float(gradient(n1-dx*direction(s)/2)@direction(s))
        restore=brentq(inverse,0.,1.,xtol=5e-15,rtol=1e-15)
        if abs(restore-settled)>2e-13 or np.max(abs(n1-dx*direction(restore)-n))>5e-15:
            raise ValueError("Reverse update did not recover the same finite event.")
        trace_flags.append(dict(extent=float(dx),frozen_weight=base,conserved_weight=settled,
                                exact_energy_residual=float(fixed),frozen_energy_drift=float(old),
                                quadratic_prediction=float(theory),reverse_weight=float(restore)))
    # A one-sided receiving bracket cannot be repaired by negative weights.
    bad=ebase.copy()
    bad[1:3]=[ed*1.04,ed*1.08]
    not_inside=all((bad+(gradient(n)-ebase))@direction(w)>0 for w in (0,1))
    if not not_inside:
        raise ValueError("One-sided rejection fixture did not exercise the guard.")
    coupled=[]
    gcap=np.array([.9,.8,.9,1.2],dtype=DT)
    is_bose=np.array([False,False,False,True])
    for tag,f in (('forward',[.30,.048,.073,.032]),('inverse',[.014,.67,.58,.31])):
        n0=gcap*np.asarray(f,dtype=DT)
        down=-DT(.8)*min(gcap[0]-n0[0],n0[1],n0[2],n0[3])
        up=DT(.8)*min(n0[0],gcap[1]-n0[1],gcap[2]-n0[2])
        def targets(dx):
            w=weight(DT(dx),n0)
            v=direction(w)
            y=n0+DT(dx)*v
            a,Hv=averaged_alphas(n0,y,gcap,is_bose)
            A0=v@a
            ad=(1-DT(w))*a[1]+DT(w)*a[2]
            s_forward=-np.logaddexp(DT(0),a[0])-np.logaddexp(DT(0),-ad)-np.log(-np.expm1(-a[3]))
            signed=np.sign(A0)*np.exp(max(s_forward,s_forward-A0))*(-np.expm1(-abs(A0)))
            return y,w,A0,signed
        for h in (.02,.08,.20):
            q=brentq(lambda d:float(d-h*targets(d)[-1]),float(down),float(up),xtol=5e-16)
            y,w,a,net=targets(q)
            err=abs(entropy(y,gcap,is_bose)-entropy(n0,gcap,is_bose)-h*a*net)
            Ee=abs(total(y)-total(n0))
            if (min(y)<=0. or max((y/gcap)[~is_bose])>=1.
                or abs(q-h*net)>3e-14 or Ee>5e-16*total(n0)
                or a*net<=0. or err>2e-14 or (q>0)!=(tag=='forward')):
                raise ValueError("Variable-weight implicit pair failed its entropy or full-energy check.")
            coupled.append(dict(type=tag,timestep=h,extent=q,finite_weight=w,
                               abs_full_energy_change=float(Ee),abs_entropy_identity_defect=float(err),
                               entropy_increase=float(h*a*net)))
    return dict(case=name,masses=[float(x),float(m),float(b)],coupling=float(g),
                manufactured_outside_bracket=bool(not_inside),controls=trace_flags,coupled_entropy= coupled)


if __name__=='__main__':
    print(json.dumps(dict(scope="Conservative finite event in a specified quadratic proxy; not a matched dispersion or physical source.",
        cases=[model('ordinary',.046,.002,.007,.36,-.33,.005),
               model('soft-scalar',.048,.0009,.003,.65,-.9998,.012)]),indent=2,allow_nan=False))
