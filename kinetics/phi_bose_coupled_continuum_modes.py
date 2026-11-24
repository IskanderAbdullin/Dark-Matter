"""Shared invariant physical χ1↔χ2+φ quantum continuum strong generators.

Modewise gain/loss for three simultaneous (possibly nonthermal) occupations.
No finite-time fully coupled integrator yet. Positivity and Pauli bounds are
inward vector-field properties; 2 numbers, energy, H follow from same events.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit
from kinetics.phi_bose_physical_cutcell_weak import M1,M2,MPHI,G1,KSTAR,ESTAR,WIDTH,shell
from kinetics.phi_bose_invariant_parent_support import parent_interval,global_support_floor

G2=4
GPHI=1


def fd1(E):return expit(-np.asarray(E)/.9+.1*np.sin(3*np.asarray(E)))
def fd2(E):return expit(-np.asarray(E)/1.+.07*np.cos(2*np.asarray(E)))
def bephi(k):return (1+.03*np.sin(3*np.asarray(k)))/np.expm1(np.asarray(k)+.005)


def parent_mode(p,f2,bphi,nk=32):
    """True Pauli gain/loss per χ1 mode at p, inverse gain G≥0 and loss L≥0."""
    E1=math.hypot(p,M1)
    low,high=shell(p)[1:3]
    if p==0:
        k=np.array([ESTAR]); weights=np.array([1.])
    else:
        z,w=leggauss(nk);ln=math.log(high/low)
        k=low*np.exp((z+1)*ln/2)
        weights=w*k*ln/(2*p*KSTAR/M1)/2
    ff2=np.asarray(f2(E1-k))
    bb=np.asarray(bphi(k))
    if np.any((ff2<0)|(ff2>1)) or np.any(bb<0):raise ValueError('invalid occupations')
    decay=WIDTH*M1/E1
    gain=decay*float(weights@(ff2*bb))
    loss=decay*float(weights@((1-ff2)*(1+bb)))
    return gain,loss


def daughter_mode(e2,P,f1,bphi,nk=32,g2=G2):
    """Pauli gain/loss per daughter χ2 energy using two-sided parent support."""
    bounds=parent_interval(e2,M2,P)
    if bounds is None or e2==M2:return 0.,0.
    pa,pb=bounds
    ka=math.hypot(pa,M1)-e2
    kb=math.hypot(pb,M1)-e2
    ka=max(ka,MPHI)
    z,w=leggauss(nk);lg=math.log(kb/ka)
    k=ka*np.exp((z+1)*lg/2)
    dk=w*k*lg/2
    f=np.asarray(f1(e2+k)); bb=np.asarray(bphi(k))
    if np.any((f<0)|(f>1)) or np.any(bb<0):raise ValueError('invalid occupations')
    p2=math.sqrt((e2-M2)*(e2+M2))
    pref=G1*WIDTH*M1*M1/(2*g2*KSTAR*e2*p2)
    return pref*float(dk@(f*(1+bb))),pref*float(dk@((1-f)*bb))


def scalar_mode(k,P,f1,f2,np_parent=32,gphi=GPHI):
    """Bose source S>=0 and signed D; dot b=S-D b with physical scalar DOS."""
    bounds=parent_interval(k,MPHI,P)
    if bounds is None:return 0.,0.
    pa,pb=bounds
    z,w=leggauss(np_parent)
    p=pa+(pb-pa)*(z+1)/2
    dp=w*(pb-pa)/2
    E1=np.hypot(p,M1)
    f=np.asarray(f1(E1));h=np.asarray(f2(E1-k))
    if np.any((f<0)|(f>1)|(h<0)|(h>1)):raise ValueError('invalid fermion occupations')
    pphi=math.sqrt((k-MPHI)*(k+MPHI))
    pref=G1*WIDTH*M1*M1/(2*gphi*KSTAR*k*pphi)
    weight=dp*p/E1
    return pref*float(weight@(f*(1-h))),pref*float(weight@(h-f))


def integrated_parent_weak(P,f1,f2,bphi,nparent=32,nangle=32):
    """Single shared positive on-shell event measure for all 3 N/Q and H.

    This yields a coupled CONTINUUM vector field identity, not a bounded
    finite-time discrete solver. The same event flux is deposited in each Q.
    """
    x,w=leggauss(nparent)
    p=P*(x+1)/2
    dp=w*P/2
    z,v=leggauss(nangle)
    ret=np.zeros(6) # Rφ,Qφ,Q1,Q2,H, min_event_weight
    ret[5]=math.inf
    for pi,dpi in zip(p,dp):
        E1=math.hypot(pi,M1)
        low,high=shell(pi)[1:3]
        lg=math.log(high/low)
        k=low*np.exp((z+1)*lg/2)
        dk=v*k*lg/2
        B=pi*KSTAR/M1
        W=G1*pi*pi*dpi/(2*math.pi**2)*WIDTH*M1/E1*dk/(2*B)
        f=float(np.asarray(f1(E1)))
        h=np.asarray(f2(E1-k));bb=np.asarray(bphi(k))
        F=f*(1-h)*(1+bb)
        R=(1-f)*h*bb
        J=F-R
        flux=W*J
        ret[0]+=float(sum(flux))
        ret[1]+=float(flux@k)
        ret[2]-=float(sum(flux))*E1
        ret[3]+=float(flux@(E1-k))
        ret[4]+=float(np.sum(flux*np.log(F/R)))
        ret[5]=min(ret[5],float(np.min(W)))
    return dict(Rphi=ret[0],Qphi=ret[1],Q1=ret[2],Q2=ret[3],
                energy_defect=sum(ret[1:4]),entropy=ret[4],minW=ret[5],
                charge_12_residual=0.,charge_1phi_residual=0.)


def integrated_scalar_weak(P,f1,f2,bphi,nk=48,nparent=36):
    """Independent scalar energy-first total Rφ and Qφ for arbitrary occupations."""
    low=global_support_floor(MPHI,P)
    high=shell(P)[2]
    cuts=[low,ESTAR,high] if low<ESTAR<high else [low,high]
    z,w=leggauss(nk)
    N=Q=0.
    for a,b in zip(cuts[:-1],cuts[1:]):
        lg=math.log(b/a)
        k=a*np.exp((z+1)*lg/2)
        dk=w*k*lg/2
        for ki,dki in zip(k,dk):
            S,D=scalar_mode(float(ki),P,f1,f2,nparent)
            bmode=float(np.asarray(bphi(float(ki))))
            source=S-D*bmode
            dos=ki*math.sqrt((ki-MPHI)*(ki+MPHI))/(2*math.pi**2)
            N+=dki*dos*source
            Q+=dki*dos*ki*source
    return N,Q


def audit():
    rows=[]
    for P in (2,8,16):
        parent=integrated_parent_weak(P,fd1,fd2,bephi,96,96)
        for n in (16,32,48,64,96):
            scalar=integrated_scalar_weak(P,fd1,fd2,bephi,n,n)
            rows.append(dict(P=P,order=n,
                             Rphi_parent=parent['Rphi'],
                             Qphi_parent=parent['Qphi'],
                             Q1_parent=parent['Q1'],
                             Q2_parent=parent['Q2'],
                             entropy_parent=parent['entropy'],
                             energy_defect_parent=parent['energy_defect'],
                             relative_Rphi_error=scalar[0]/parent['Rphi']-1,
                             relative_Qphi_error=scalar[1]/parent['Qphi']-1))
    return dict(schema_version=1,
                scope="arbitrary nonthermal 3-species CONTINUUM mode generators and weak event ledgers only; no full coupled finite-time integrator",
                rows=rows,class_ready=False)

if __name__=='__main__':
    import json
    print(json.dumps(audit(),indent=2,allow_nan=False))
