"""Projected k^2 attenuation of the locally locked static acoustic pair."""
import argparse
import json
from pathlib import Path

import numpy as np
from charged_lepton_portal_local_longitudinal_reference import model
from positive_fd_portal_quadrature_reference import PositivePortal


def analyze(T,xi,portal,r,scenario):
    L,H,*_=model(T,xi,portal,r,0.,scenario)
    L1=model(T,xi,portal,r,1.,scenario)[0]
    hs,V=np.linalg.eigh(H)
    W=(V*np.sqrt(hs))@V.T
    J=(V/np.sqrt(hs))@V.T
    A=(W@L@J).real
    B=(-(W@(L1-L)@J).imag)
    if max(np.max(abs(A-A.T)),np.max(abs(B-B.T)))>3e-14:
        raise ValueError("Thermodynamic whitened decay/streaming matrices must be symmetric.")
    rates,U=np.linalg.eigh(-.5*(A+A.T))
    if np.count_nonzero(rates>1e-11)!=3 or min(rates)<-1e-13:
        raise ValueError("Expected two conserved and three dissipative thermodynamic modes.")
    N,F=U[:,:2],U[:,2:]
    slow=N.T@B@N
    speeds,P=np.linalg.eigh(slow)
    v=N@P[:,1]
    if not (speeds[1]>0. and abs(speeds.sum())<2e-13):
        raise ValueError("No equal and opposite acoustic pair.")
    g=B@v
    comp=(F.T@g)/np.sqrt(rates[2:])
    attenuation=float(comp@comp)
    gs=g.copy(); gs[3:]=0.
    gv=g-gs
    cs=(F.T@gs)/np.sqrt(rates[2:])
    cv=(F.T@gv)/np.sqrt(rates[2:])
    if abs(cs@cv)>2e-11*max(attenuation,1e-20):
        raise ValueError("Scalar and vector decay did not separate by parity.")
    if np.linalg.norm(F.T@gv)>2e-12*max(1.,np.linalg.norm(B)):
        raise ValueError("A locally equilibrated compression drove relative velocity at first order.")
    from_scalar=F@((F.T@gs)/rates[2:])
    through_vector=F@(F.T@(B@from_scalar))
    through_vector[:3]=0.
    q4_numerator=float(through_vector@through_vector)
    limits=[]
    for portion in (.09,.045,.0225):
        k=portion*min(rates[2:])
        M=model(T,xi,portal,r,k,scenario)[0]
        w,z=np.linalg.eig(W@M@J)
        idx=int(np.argmin(abs(w+1j*k*speeds[1])))
        q=z[:,idx]
        f=F.T@q
        inferred=float(np.sum(abs(f)**2*rates[2:])/np.vdot(q,q).real)
        if abs(inferred+w[idx].real)>3e-13:
            raise ValueError("Entropy quotient disagrees with direct real eigenvalue.")
        limits.append(dict(k_over_heat_rate=k,damping_over_k_squared=inferred/(k*k),phase_speed=float(-w[idx].imag/k)))
    if abs(limits[-1]["damping_over_k_squared"]/attenuation-1)>2e-3:
        raise ValueError("Exact acoustic mode did not approach projected damping.")
    return dict(scenario=scenario or "full_momentum",sound_speed=float(speeds[1]),
                attenuation_coefficient=attenuation,scalar_decay_term=float(cs@cs),vector_decay_term=float(cv@cv),
                decay_rates=rates[2:].tolist(),vector_return_numerator=q4_numerator,
                spectral_small_k_checks=limits)


def run(receipt):
    p=PositivePortal((13,36,26),tail=65.)
    out=[]
    for r in receipt['rows']:
        for xi in (1.,10.):
            modes=[analyze(r['T'],xi,p,r,s) for s in (None,"heat_as_momentum","ignore_elastic_momentum")]
            exact=modes[0]
            if any(abs(d['sound_speed']-exact['sound_speed'])>2e-13 or
                   abs(d['scalar_decay_term']-exact['scalar_decay_term'])>2e-13 for d in modes[1:]):
                raise ValueError("Changing pure vector decay altered locked thermodynamics or scalar loss.")
            verify=[]
            if r['T'] in (.015,.05):
                T=r['T']
                models=(None,"heat_as_momentum")
                rate=[model(T,xi,p,r,0.,s)[4] for s in models]
                coefficient=exact["vector_return_numerator"]*(1/rate[0]-1/rate[1])
                for k in (.03,.015,.0075):
                    vals=[]
                    for s in models:
                        w=model(T,xi,p,r,k,s)[2]
                        vals.append(float(w[int(np.argmax(w.imag))].real))
                    observe=(vals[0]-vals[1])/k**4
                    verify.append(dict(k_over_heat_rate=k,projected_k4_difference=coefficient,
                                       direct_difference_over_k4=observe))
                if abs(verify[-1]["direct_difference_over_k4"]/coefficient-1)>1e-2:
                    raise ValueError("Drag replacement did not enter the acoustic attenuation at fourth order.")
            out.append(dict(T=r['T'],visible_to_dark_enthalpy_ratio=xi,controls=modes,
                            fourth_order_drag_difference=verify))
    return dict(scope="Small-k acoustic expansion of the static entropy-checked thermal moment system only.",
                rows=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt')
    a=p.parse_args()
    print(json.dumps(run(json.loads(Path(a.receipt).read_text())),indent=2,allow_nan=False))
