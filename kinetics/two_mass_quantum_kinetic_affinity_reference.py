"""Stable signed Bose/Pauli event in the near-rest kinetic coordinate."""
import json
import numpy as np
from scipy.special import roots_legendre
from two_mass_rest_strip_exact_area_reference import sinhm
from two_mass_kinetic_offset_deposition_reference import event, DT


def evaluate(ks,kf,Ts,Tf,Tp,etas,etaf,rest_affinity):
    ks,kf,Ts,Tf,Tp,etas,etaf,rest_affinity=map(DT,(ks,kf,Ts,Tf,Tp,etas,etaf,rest_affinity))
    if min(Ts,Tf,Tp)<=0 or min(ks,kf)<0:
        raise ValueError("Use positive thermal scales and nonnegative kinetic offsets.")
    a=ks/Ts-etas
    b=kf/Tf-etaf
    c=(ks+kf)/Tp-etas-etaf-rest_affinity
    if a<=0:
        raise ValueError("Nonpositive Bose exponent.")
    A=rest_affinity+ks*(Tp-Ts)/(Tp*Ts)+kf*(Tp-Tf)/(Tp*Tf)
    fwd=-np.logaddexp(DT(0),c)-np.logaddexp(DT(0),-b)-np.log(-np.expm1(-a))
    rev=fwd-A
    S=np.sign(A)*np.exp(max(fwd,rev))*(-np.expm1(-abs(A)))
    if S*A<0 or (S==0)!=(A==0):
        raise ValueError("Sign or exact equilibrium failed.")
    return float(S),float(A),float(fwd)


def reference(v):
    import mpmath as mp
    with mp.workdps(100):
        s,f,Ts,Tf,Tp,u,z,r=map(lambda y:mp.mpf(float(y)),v)
        f_B=1/mp.expm1(s/Ts-u)
        f_F=1/(1+mp.exp(f/Tf-z))
        f_P=1/(1+mp.exp((s+f)/Tp-u-z-r))
        return float(f_P*(1-f_F)*(1+f_B)-(1-f_P)*f_F*f_B)


def event_checks():
    _,K,_=event(2e-23,.061,.122,4e-12,-.31)
    rows=[]
    for Ts,Tf,Tp,x in ((.03,.03,.03,0.),(.03,.03,.03,1e-28),(.03,.03,.03,-1e-28),
                       (.03*(1+1e-6),.03,.03,0.),(.03,.03*(1-2e-6),.03,0.)):
        v=(*map(float,K),Ts,Tf,Tp,-2e-23,-.35,x)
        S,A,F=evaluate(*v)
        known=reference(v)
        error=abs(S-known) if A==0 else abs(S-known)/abs(S)
        if error>(1e-75 if A==0 else 2e-13):
            raise ValueError(f"Independent reference failed: {v}, {S}, {known}")
        rows.append(dict(kinetic_inputs=v,event_affinity=A,stable_net=S,mp_net=known,
                         conventional_float_subtract=float(np.exp(np.float64(F))-np.exp(np.float64(F-A))),
                         reference_error=error,reference_error_absolute_when_equilibrated=(A==0)))
    return rows


def strip_check(kind, edge, order=32):
    # e, boson mass, fermion mass and thermal parameters are manufactured.
    e,m,b,T=map(DT,(.031,.065,.018,.03))
    Y=2*np.arcsinh(np.sqrt(e*(2*(m+b)+e)/(4*m*b)))
    D=m*b*np.sinh(Y)
    z,w=(np.array(v,dtype=DT) for v in roots_legendre(order))
    p=edge*m*(1+z)/2
    dp=edge*m*w/2
    if kind=='affinity_tuned':
        Ts,Tf,Tp=T,T,T*DT(1.03)
        # Only the boson temperature differs from the parent, so affinity
        # vanishes proportionally to scalar kinetic excess at p=0.
        Tf=Tp
        eta,rest=DT(0),DT(0)
    else:
        Ts=Tf=Tp=T
        eta=DT(-.002) if kind=='positive_bose_exponent_floor' else DT(0)
        rest=DT(1e-7)
    minor=DT(-.35)
    answer=DT(0)
    for q,h in zip(p,dp,strict=True):
        u=np.arcsinh(q/m)
        ks=q*q/(np.hypot(q,m)+m)
        other_center=b*(2*np.sinh(Y/2)**2+2*np.cosh(Y)*np.sinh(u/2)**2)
        width=2*b*np.sinh(Y)*np.sinh(u)
        line=DT(0)
        for zz,ww in zip(z,w,strict=True):
            kf=other_center+width*zz/2
            S,_,_=evaluate(ks,kf,Ts,Tf,Tp,eta,minor,rest)
            line+=DT(S)*ww/2
        answer+=line*width*h*q/np.hypot(q,m)
    cut=edge*m
    area=D*sinhm(2*np.arcsinh(edge))/2
    if kind=='positive_bose_exponent_floor':
        upper=area*min(1,abs(rest))/(-np.expm1(eta))
    elif kind=='affinity_tuned':
        upper=area*(cut*cut/(np.hypot(cut,m)+m)+Ts)*abs(Tp-Ts)/(Tp*Ts)
    else:
        # Integral of (1 + T/K_s) over dE_s dE_f.
        upper=min(1,abs(rest))*(area+Ts*2*D/m**2*(cut+m*np.arcsinh(edge)))
    if not 0<answer<=upper:
        raise ValueError("Thermal strip source escaped its independent positive bound.")
    return dict(p_over_scalar_mass=float(edge),integrated_signed_source=float(answer),
                independent_bound=float(upper),geometric_area=float(area))


def run():
    profiles=[]
    for name,power in (('saturated_nonzero_affinity',1),
                       ('positive_bose_exponent_floor',3),('affinity_tuned',3)):
        rows=[]
        for p in (1e-2,1e-4,1e-5,1e-6,1e-8,1e-10):
            low,high=[strip_check(name,DT(p),n) for n in (16,32)]
            a='integrated_signed_source'
            if abs(low[a]/high[a]-1)>3e-10:
                raise ValueError("Independent order check failed.")
            rows.append(high)
        z=[r['integrated_signed_source'] for r in rows]
        slope=np.log(z[-1]/z[-2])/np.log(rows[-1]['p_over_scalar_mass']/rows[-2]['p_over_scalar_mass'])
        if abs(slope-power)>2e-5:
            raise ValueError("Near-rest power was not resolved.")
        profiles.append(dict(profile=name,near_rest_power=float(slope),rows=rows))
    return dict(scope="Thermal signed Bose/Pauli event and bounded omitted strip, not physical rates.",
                events=event_checks(),profiles=profiles)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
