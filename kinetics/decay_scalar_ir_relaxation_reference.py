"""Momentum-resolved IR clock for chi1 <-> chi2 + phi in thermal fermion baths."""
import argparse
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit
from three_species_decay_onsager_reference import M,m,mp,P,Es,width,susceptibility
from sm_dark_entropy_exchange_reference import StandardModel,MPL

def limits(p):
    assert p>=0.
    e=np.hypot(p,mp)
    lo=M*(Es*Es+p*p)/(e*Es+p*P)
    de=2*M*p*P/(mp*mp)
    return e,lo,de

def log_rate(p,T,fermion_eta=0.):
    e,L,ell=limits(p)
    A=(L-e)/T-fermion_eta
    b=e/T
    span=ell/T
    # Divide by the interval before taking p->0: the exact kinematic
    # width cancels the apparent 1/p in the scalar-bin weak rate.
    if span<1e-4:
        j=quad(lambda u:np.exp(-span*u)*expit(A+span*u)*expit(A+b+span*u),
               0.,1.,epsabs=1e-13,epsrel=1e-13)[0]
        norm=np.log(j)
    else:
        j=quad(lambda u:np.exp(-u)*expit(A+u)*expit(A+b+u),0.,min(55.,span),
               epsabs=1e-13,epsrel=1e-13)[0]
        norm=np.log(j)-np.log(span)
    return np.log(4*width)+3*np.log(M)-2*np.log(mp)-np.log(e)-A+np.log(-np.expm1(-b))+norm

def validate_all(T,p):
    y=log_rate(p,T)
    e,L,ell=limits(p)
    def evaluate(fphi):
        # Direct weak source from the two separate occupation products
        # (close to each other at equilibrium). Ignore negligible tails.
        def source(z):
            f1=expit(-(L+T*z)/T)
            f2=expit(-(L-e+T*z)/T)
            return f1*(1-f2)*(1+fphi)-f2*(1-f1)*fphi
        v=quad(source,0.,min(55.,ell/T),epsabs=1e-150,epsrel=5e-10)[0]*T
        return 4*width*M*M/(2*P*p*e)*v
    f=1/np.expm1(e/T)
    d=.02*f
    slope=(evaluate(f-d)-evaluate(f+d))/(2*d)
    assert slope>0
    return abs(y-np.log(slope))

def main(table):
    eos=StandardModel(table)
    rows=[]
    for T in (.01,.02,.1,.3,.6882487043,.95):
        ds=sum(susceptibility(T,mi,g,BE)[1] for mi,g,BE in ((M,4,False),(m,4,False),(mp,1,True)))
        H=np.sqrt((eos.state(T)[0]+ds)/3)/MPL
        shift=np.log(H)
        # Search ascending momentum for the first local damping/H crossing.
        t=np.geomspace(1e-8,1.,160)
        r=np.array([log_rate(x*T,T)-shift for x in t])
        roots=[np.exp(brentq(lambda l:log_rate(np.exp(l)*T,T)-shift,np.log(x),np.log(z)))
                for x,z,a,b in zip(t[:-1],t[1:],r[:-1],r[1:]) if a<0<=b]
        assert len(roots)==1
        r0=log_rate(0.,T)-shift
        rnear=log_rate(1e-16*mp,T)-shift
        assert abs(rnear-r0)<1e-8
        eps=max(validate_all(T,q*T) for q in (.2,.5,1.))
        assert eps<1e-8
        qmin=roots[0]
        rph=limits(0.)[1]/T
        rows.append(dict(T_GeV=T,log10_rest_mode_damping_over_H=r0/np.log(10.),
                          first_crossing_p_over_T=qmin,first_crossing_p_GeV=qmin*T,
                          parent_threshold_at_rest_over_T=rph,
                          local_detailed_balance_slope_log_error=eps,
                          example_log10_rates_over_H={str(q):(log_rate(q*T,T)-shift)/np.log(10.) for q in (1e-6,.0001,.001,.01,.1,1.)}))
    return dict(rows=rows,scope="Decay/inverse decay channel only, thermal fermion baths. An integrated equilibration clock cannot ensure the low-momentum scalar tail equilibrates. No condensate or total freeze-out claim.")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("eos")
    print(json.dumps(main(p.parse_args().eos),indent=2,allow_nan=False))
