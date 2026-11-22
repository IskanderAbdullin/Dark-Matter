"""Exact-energy elementary χ1<->χ2+real φ Pauli/Bose event, bounded implicit.

Scoped conservative quantum reaction test, not a continuum physical Γ/time
operator or an abundance/IR/FRW/86eV model. All energies must be verified
against a Lorentz-allowed moving on-shell decay outside this module.
"""
from __future__ import annotations
import math
from scipy.optimize import brentq
from scipy.special import xlogy


def fermion_entropy(f):
    if not 0<=f<=1:raise ValueError('Fermi bounds')
    return -float(xlogy(f,f)+xlogy(1-f,1-f))


def bose_entropy(b):
    if b<0:raise ValueError('negative Bose mode')
    return float(xlogy(1+b,1+b)-xlogy(b,b))


def quantum_J(f1,f2,b):
    return f1*(1-f2)*(1+b)-(1-f1)*f2*b


def entropy(f1,f2,b,w1,w2,wphi):
    return w1*fermion_entropy(f1)+w2*fermion_entropy(f2)+wphi*bose_entropy(b)


def ledger(f1,f2,b,E1,E2,Ephi,w1,w2,wphi):
    return {'Nfermion':w1*f1+w2*f2,
            'N1_plus_phi':w1*f1+wphi*b,
            'energy':w1*f1*E1+w2*f2*E2+wphi*b*Ephi,
            'entropy':entropy(f1,f2,b,w1,w2,wphi)}


def implicit_event(f1,f2,b,*,E1,E2,Ephi,w1,w2,wphi,W,dt):
    vals=(f1,f2,b,E1,E2,Ephi,w1,w2,wphi,W,dt)
    if not all(math.isfinite(x) for x in vals) or not(0<=f1<=1 and 0<=f2<=1 and b>=0 and min(w1,w2,wphi,W)>0 and dt>=0):
        raise ValueError('invalid Pauli/Bose reaction')
    if abs(E1-E2-Ephi)>2e-12*max(E1,E2,Ephi,1):
        raise ValueError('off shell: E1 must equal E2+Ephi')
    if dt==0:return f1,f2,b,0.
    lower=max((f1-1)*w1,-f2*w2,-b*wphi)
    upper=min(f1*w1,(1-f2)*w2)
    def flux(x):
        return quantum_J(f1-x/w1,f2+x/w2,b+x/wphi)
    def fun(x):return x-dt*W*flux(x)
    # The physical endpoints bracket the unique quantum reaction extent.
    x=brentq(fun,lower,upper,xtol=1e-18,rtol=9e-16)
    after=(f1-x/w1,f2+x/w2,b+x/wphi)
    old=ledger(f1,f2,b,E1,E2,Ephi,w1,w2,wphi)
    new=ledger(*after,E1,E2,Ephi,w1,w2,wphi)
    if new['entropy']<old['entropy']-1e-13:
        raise ArithmeticError('entropy decrease')
    return *after,x
