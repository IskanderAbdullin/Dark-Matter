"""Conditional anomaly-free L_e−L_mu vector portal lepton-channel heat audit.

Gauged U(1)_{L_e−L_mu} is an explicit electroweak-covariant UV *candidate*
(not a demonstrated embedding of χ1↔χ2+φ). Dirac vector-like χ couples
to the new gauge boson with charge qχ; eL,eR,νeL +g', μL,μR,νμL
−g'. For a mediator mass large compared with all t and T:
Gχℓ≈gχ qχ g'/mV², same |G| for electrons and both ν species.
Physical ν are one helicity each, plus antineutrinos. A vector χ
symmetric trace kills antisymmetric neutrino axial interference,
leaving same spin-averaged massless vector kernel per helicity.
All heat values are CONDITIONAL 5D tree-contact FD equal-Te/Tν.
No real massive mediator denominator, bounds, ν decoupling, cosmology.
"""
from __future__ import annotations

import argparse
import json
import math
import numpy as np

from kinetics.independent_massive_chie_5d import cm_heat


def gamma_matrices():
    I=np.eye(2,dtype=complex)
    z=np.zeros((2,2),complex)
    p=[np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),
       np.array([[1,0],[0,-1]],complex)]
    g=[np.block([[I,z],[z,-I]])]
    g += [np.block([[z,s],[-s,z]]) for s in p]
    g5=1j*g[0]@g[1]@g[2]@g[3]
    return g,(np.eye(4)-g5)/2


def massive_chi_neutrino_trace(M,p,c):
    """Compare explicit ν_L one-helicity vs Dirac-vector unpolarized trace."""
    g,PL=gamma_matrices()
    eta=[1,-1,-1,-1]
    one=np.eye(4,dtype=complex)
    def slash(k):
        return sum(eta[i]*k[i]*g[i] for i in range(4))
    E=math.hypot(p,M)
    si=math.sqrt(1-c*c)
    a=np.array([E,0,0,p])
    b=np.array([p,0,0,-p])
    ap=np.array([E,p*si,0,p*c])
    bp=np.array([p,-p*si,0,-p*c])
    chi=[[np.trace((slash(ap)+M*one)@g[u]@
                    (slash(a)+M*one)@g[v])/2
          for v in range(4)] for u in range(4)]
    vector=[[np.trace(slash(bp)@g[u]@slash(b)@g[v])/2
             for v in range(4)] for u in range(4)]
    left=[[np.trace(slash(bp)@g[u]@PL@slash(b)@g[v]@PL)
           for v in range(4)] for u in range(4)]
    def contract(mat):
        return sum(eta[i]*eta[j]*chi[i][j]*mat[i][j]
                   for i in range(4) for j in range(4))
    v=contract(vector)
    l=contract(left)
    s=M*M+2*p*(E+p)
    t=-2*p*p*(1-c)
    u=2*M*M-s-t
    expected=2*((s-M*M)**2+(u-M*M)**2+2*M*M*t)
    return float(v.real),float(l.real),expected,float(abs(v.imag)+abs(l.imag))


def heat_lepton_comparison(T_me,*,M_me=.1,G_unit=.001,Tchi_over_T=.9,
                           n=28,cut_kinetic_over_T=30):
    """Independent massive vs massless 5D Pauli CM heat, SAME microphysics.

    Units me=1; Q in me^5 with dimensionless G=.001 me^-2.
    Neutrino gν=4 means νe, νbar_e, νμ, νbar_μ (one helicity each).
    Electron ge=4 means e−, e+ with two helicities each.
    """
    if not (T_me>0 and M_me>0 and G_unit>0 and Tchi_over_T>0
            and isinstance(n,int) and n>=12 and cut_kinetic_over_T>0):
        raise ValueError("positive benchmark scales")
    T=T_me
    Tc=Tchi_over_T*T
    M=M_me
    pmax=math.sqrt((M+cut_kinetic_over_T*T)**2-M*M)
    vals={}
    for mass,label in ((0.,"two_neutrino_flavors"),(1.,"electron_positron")):
        kmax=math.sqrt((mass+cut_kinetic_over_T*T)**2-mass*mass)
        vals[label]=cm_heat(M,mass,G_unit,T,Tc,
                            pmax=pmax,kmax=kmax,
                            nr=n,nang=n,nphi=n,
                            gchi=4.,ge=4.)
    return {"T_over_me":T_me,"Mchi_over_me":M_me,
            "Tchi_over_Tlepton":Tchi_over_T,
            "G_me_minus2":G_unit,
            "Q_chi_two_active_nu_G2_me9_scaled":vals["two_neutrino_flavors"]/G_unit**2,
            "Q_chi_e_minus_e_plus_G2_me9_scaled":vals["electron_positron"]/G_unit**2,
            "nu_over_e_heat_ratio":vals["two_neutrino_flavors"]/vals["electron_positron"] if vals["electron_positron"] else None,
            "g_active_nu":4,"g_electron_positron":4,
            "assumed_nu_equal_e_temperature":True,
            "no_physical_gamma_or_86eV_relic":True}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--T-over-me",type=float,default=.1)
    ap.add_argument("--n",type=int,default=28)
    op=ap.parse_args()
    print(json.dumps(heat_lepton_comparison(op.T_over_me,n=op.n),indent=2))
