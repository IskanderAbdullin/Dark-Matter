"""Conformal-time increments on the explicitly modeled cooling intervals."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import constants
from scipy.integrate import quad

from archived_collisionless_streaming_reference import Background
from finite_cooling_reheat_portal_reference import TF, SourceFold


def run(eos,reference):
    ref=json.loads(Path(reference).read_text())
    hexdigest=hashlib.sha256(Path(eos).read_bytes()).hexdigest()
    assert ref["eos_sha256"]==hexdigest
    B=Background(eos)
    K=SourceFold(eos,(4,9,8),65,4)
    ai=float(B.scale(TF))
    z=np.log(TF)
    units=constants.hbar*constants.c/(1e9*constants.e*1e6*constants.parsec)
    kmax=10*B.H0/100
    output=[]
    for v in (ref["controls"][-1]["rows"][0],ref["controls"][-1]["rows"][-1]):
        tr,top=v["T_R_GeV"],v["Tmax_GeV"]
        r=np.log(tr)
        ar=float(B.scale(tr))
        ratio=top/tr
        aratio=ratio**(-8/3)
        pre=2*(-np.expm1(.5*np.log(aratio)))*units/(ar*K.H(r))
        def pre_integrand(v):
            return units*np.exp(.5*v)/(ar*K.H(r))
        pre_direct=quad(pre_integrand,np.log(aratio),0.,epsabs=1e-19,epsrel=1e-12)[0]
        assert abs(pre-pre_direct)<1e-17
        def radiation(u):
            a=ai*(K.s(z)/K.s(u))**(1/3)
            return (1+K.grad(u)/3)*units/(a*K.H(u))
        knots=K.bg.sg.x[(K.bg.sg.x>z)&(K.bg.sg.x<r)]
        post=quad(radiation,z,r,points=knots,epsabs=1e-19,epsrel=1e-10,limit=2000)[0]
        # Use separate scale-to-photon-T inverse to test the log T-to-a measure.
        from scipy.optimize import brentq
        def in_lna(v):
            u=brentq(lambda t: 3*(np.log(ai)-v)+np.log(K.s(z)/K.s(t)),z,r,xtol=5e-14)
            return units/(np.exp(v)*K.H(u))
        val=quad(in_lna,np.log(ar),np.log(ai),epsabs=1e-15,epsrel=1e-9,limit=300)[0]
        assert abs(post/val-1)<3e-8
        output.append(dict(T_R_GeV=tr,Tmax_GeV=top,a_at_TR=ar,a_at_1MeV=ai,a_at_Tmax=ar*aratio,
            Delta_tau_prescribed_pre_cooling_Mpc=float(pre),Delta_tau_prescribed_radiation_Mpc=float(post),
            kmax_times_Delta_tau_on_modeled_intervals=float(kmax*(pre+post)),
            log_scale_clock_relative_error=abs(post/val-1.)))
    return dict(eos_sha256=hexdigest,kmax_h_per_Mpc=10.,cases=output,
            scope="Conformal intervals only, from the prescribed Tmax start to 1 MeV. Not the absolute conformal time from the big bang or a construction of perturbations in the omitted heating era. The main tracer uses an independent adiabatic CLASS reference after 1 MeV.")


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("reference")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference),indent=2,allow_nan=False))
