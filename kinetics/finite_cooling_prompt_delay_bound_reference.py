"""Vacuum heavy-parent delay bound for the finite cooling reference."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from conditional_relic_normalization_reference import conversion
from finite_cooling_reheat_portal_reference import TF, SourceFold
from three_species_decay_onsager_reference import M, P, mp, width


def phase_integral(K,tr,tmax):
    r,h=np.log([tr,tmax])
    L=np.log(TF)
    out=[]
    for i,(a,b) in enumerate(((L,r),(r,h))):
        v=np.zeros(3)
        if b<=a:
            out.append(v)
            continue
        def func(t):
            f=np.empty_like(t)
            if i==0:
                f[:]=(3+K.grad(t))/(3*K.H(t)*K.s(t))
            else:
                f[:]=(8/3)*np.exp(-12*(t-r))/(K.H(r)*K.s(r))
            v=[]
            for u,w in zip(t.flat,f.flat):
                E,G=K.C[1].particles(np.exp(u),6)
                b=(E/M)**2
                v.append(w*np.array([G.sum(),G@b,G@(b*b)]))
            return np.asarray(v).reshape(*t.shape,3)
        out.append(K.fold(a,b,np.r_[np.linspace(a,b,39),K.bg.sg.x],func))
    return np.array(out)


def run(eos,reference):
    ref=json.loads(Path(reference).read_text())
    eos_sha=hashlib.sha256(Path(eos).read_bytes()).hexdigest()
    assert ref["eos_sha256"]==eos_sha
    rows=ref["controls"][-1]["rows"]
    K=SourceFold(eos,(8,23,17),115,6)
    K2=SourceFold(eos,(6,18,14),85,4)
    L=np.log(TF)
    factor=conversion(eos)["Omega_h2_per_GeV_times_yield"]
    out=[]
    for row in (rows[0],rows[-1]):
        tr,tmax=row["T_R_GeV"],row["Tmax_GeV"]
        r=np.log(tr)
        Hmax=K.H(r)*(tmax/tr)**4
        # In the radiation segment s must increase with T (positive heat
        # capacity), and rho must increase with T. Check full monotone PCHIP
        # derivative minima panel by panel as a sign check.
        for f,k in ((K.bg.rg,4.),(K.bg.sg,3.)):
            x=np.unique(np.clip(np.r_[f.x,f.derivative(2).roots(extrapolate=False),L,r],L,r))
            assert np.min(f.derivative()(x[np.isfinite(x)])+k)>0.
        V=phase_integral(K,tr,tmax).sum(axis=0)
        cmp=phase_integral(K2,tr,tmax).sum(axis=0)
        assert np.max(abs(cmp/V-1.))<3e-5
        count_defect=abs(V[0]*factor*mp/row["scalar_prompt_fixed_count_Omega_h2"]-1.)
        assert count_defect<6e-5
        dat=[]
        for delta in (1e-9,1e-8,1e-7):
            threshold=brentq(lambda u: np.log(K.s(u)/K.s(L))-3*delta,L,min(r,L+1e-5),xtol=1e-15)
            late=phase_integral(K,tr,tr,lo=L,nt=6)[0,0] if threshold==r else                  K.fold(L,threshold,K.bg.sg.x,lambda ts: np.array([
                    [K.C[1].particles(np.exp(u),4)[1].sum()*(3+K.grad(u))/(3*K.H(u)*K.s(u))]
                      for u in ts.flat]).reshape(*ts.shape,1))[0]
            A=width/Hmax*delta
            bounds=V[1:]/V[0]*(np.array([2.,4.])/(np.e*A))**np.array([2,4])
            dat.append(dict(delta_ln_a=delta,
                fraction_born_after_one_delay_before_anchor=float(late/V[0]),
                fraction_earlier_parents_delayed_longer_than_delta_bounds_n2_n4=bounds.tolist(),
                missing_decay_fraction_upper_at_anchor=float(min(1.,late/V[0]+min(bounds))),
                absolute_p_over_T_difference_upper_for_prompt_decay_class=float(P/TF*np.expm1(delta))))
        out.append(dict(T_R_GeV=tr,Tmax_GeV=tmax,H_max_GeV=Hmax,width_over_H_max=width/Hmax,
                        parent_E_over_M_moments_2_4=(V[1:]/V[0]).tolist(),
                        relative_prompt_scalar_count_defect=float(count_defect),
                        max_relative_quadrature_change=float(np.max(abs(cmp/V-1.))),bounds=dat))
    return dict(width_GeV=width,eos_sha256=eos_sha,scope="Vacuum exponential decay of freely redshifting heavy parents on the prescribed cooling clock. Weighted moment bounds refer to the integrated quadrature source, not an independently bounded physical reheating solution.",cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    p.add_argument("reference")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.reference),indent=2,allow_nan=False))
