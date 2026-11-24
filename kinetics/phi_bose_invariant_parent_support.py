"""General two-sided invariant inverse parent support for massive 1→2 daughters.

Pure Lorentz kinematics, independent of FD/BE occupation, collision rates or
a condensate ansatz. Both fermion daughter χ2 and ultra-light scalar φ can
have FINITE two-sided allowed parent momentum intervals.
"""
import math
import json
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,MPHI,KSTAR,shell,
)


def parent_interval(daughter_energy,daughter_mass,parent_pmax):
    """Return actual physical [parent p_low,p_high] for a massive daughter.

    Daughter rest-frame momentum q*=KSTAR; define ξ=asinh(q*/m),
    lab daughter ζ=asinh(pd/m), parent y=asinh(p/M1).
    Boosted energy shell: |ξ−y|<=ζ<=ξ+y, i.e.
    |ξ−ζ|<=y<=ξ+ζ, further intersect 0<=y<=asinh(P/M1).
    In contrast to the lower-threshold-only scheme, high-parent support
    can end below P even for daughter E above its rest mass.
    """
    if daughter_mass<=0 or parent_pmax<=0:
        raise ValueError("Positive mass and parent momentum cutoff required")
    if daughter_energy<daughter_mass:return None
    q=math.sqrt(max(0.,(daughter_energy-daughter_mass)*
                          (daughter_energy+daughter_mass)))
    zeta=math.asinh(q/daughter_mass)
    xi=math.asinh(KSTAR/daughter_mass)
    ymax=math.asinh(parent_pmax/M1)
    ylo=abs(xi-zeta)
    if ylo>ymax:
        if ylo-ymax>2e-13:return None
        ylo=ymax  # floating roundoff at exact finite-parent support edge
    yhi=min(ymax,xi+zeta)
    return M1*math.sinh(ylo),M1*math.sinh(yhi)


def global_support_floor(daughter_mass,pmax):
    """Actual daughter energy min for ALL parent p∈[0,P].

    p*=M1*q*/m gives a daughter AT REST. If P<p* minimum is
    the anti-parallel energy at P. If P>=p*, minimum is its own rest mass,
    NOT the anti-parallel endpoint value at P (which rises again).
    """
    if daughter_mass<=0 or pmax<=0:
        raise ValueError("Positive mass/cutoff required")
    xi=math.asinh(KSTAR/daughter_mass)
    yp=math.asinh(pmax/M1)
    return daughter_mass*math.cosh(max(0.,xi-yp))


def audit():
    pstar_phi=M1*KSTAR/MPHI
    return dict(scalar_rest_parent_GeV=pstar_phi,
                scalar_global_min_P16_GeV=global_support_floor(MPHI,16.),
                phi_86p1eV_parent_interval_P60TeV_GeV=
                   parent_interval(86.1e-9,MPHI,60000.),
                daughter2_rest_parent_GeV=M1*KSTAR/M2,
                daughter2_Estar_parent_interval_P2_GeV=
                   parent_interval((M1*M1+M2*M2-MPHI*MPHI)/(2*M1),
                                   M2,2.),
                unsupported_86eV_rest_at_P16=True)


if __name__=="__main__":
    print(json.dumps(audit(),indent=2))
