"""Positive cut-cell Yukawa inverse/forward Bose-Pauli weak collision reference.

Integrates exactly the scalar-energy interval intersecting each parent
boosted on-shell shell, instead of asking sparse fixed Gauss angles to
hit a narrow daughter cell. Fixed thermal occupations; no strong evolution.
"""
import math
from scipy.integrate import quad
from scipy.special import expit
from kinetics.phi_bose_ir_accessibility import (
    M1, M2, MPHI, rest_kinematics, parent_threshold_for_daughter_energy,
)
from kinetics.phi_bose_sparse_angle_false_zero import first_cell_edge

Y = .001
G1 = 4
T1 = .9
T2 = 1.
TPHI = 1.


def width():
    lam=(M1*M1-(M2+MPHI)**2)*(M1*M1-(M2-MPHI)**2)
    return Y*Y*((M1+M2)**2-MPHI**2)*math.sqrt(lam)/(16*math.pi*M1**3)


def integrate(pmax=2., scalar_max=None, energy_weight=False):
    """Signed net scalar source or physical scalar heat in GeV^4 / GeV^5.

    scalar_max=None: all daughter energies; otherwise integrate only
    mphi<=Ephi<=scalar_max. Two nested adaptive QUAD integrals.
    """
    if pmax<=0: raise ValueError("positive parent cutoff")
    if scalar_max is not None and scalar_max<=MPHI: return 0.
    estar,qstar=rest_kinematics()
    pstart=0. if scalar_max is None else parent_threshold_for_daughter_energy(
        scalar_max)
    if pstart>=pmax:return 0.
    gam=width()*G1/(2*math.pi**2)
    def parent(p):
        if p==0: return 0.
        e=math.hypot(p,M1)
        b=p*qstar/M1
        kmin=(MPHI*MPHI*p*p+M1*M1*estar*estar)/(M1*(e*estar+p*qstar))
        kmax=(e*estar+p*qstar)/M1
        khi=kmax if scalar_max is None else min(kmax,scalar_max)
        if khi<=kmin:return 0.
        f1=expit(-e/T1)
        def flux(k):
            f2=expit(-(e-k)/T2)
            bose=1/math.expm1(k/TPHI)
            net=f1*(1-f2)*(1+bose)-(1-f1)*f2*bose
            return k*net if energy_weight else net
        integral=quad(flux,kmin,khi,epsabs=1e-13,
                      epsrel=1e-10,limit=200)[0]/(2*b)
        return gam*p*p*M1/e*integral
    return quad(parent,pstart,pmax,epsabs=1e-22,
                epsrel=1e-9,limit=300)[0]


def audit():
    p=2.
    top=first_cell_edge(p,16)
    cell=integrate(p,top)
    heat=integrate(p,top,True)
    total=integrate(p)
    return {"parent_pmax_GeV":p,"scalar_N":16,
            "first_cell_top_GeV":top,
            "first_cell_parent_threshold_GeV":parent_threshold_for_daughter_energy(top),
            "first_cell_Rphi_GeV4":cell,
            "first_cell_Qphi_GeV5":heat,
            "total_Rphi_GeV4":total,
            "first_cell_fraction_of_total":cell/total,
            "previous_12x12_Rphi_GeV4":-1.1649301286296807e-11,
            "old_12x12_fractional_total_underestimate":
                (total-(-1.1649301286296807e-11))/total}


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2))
