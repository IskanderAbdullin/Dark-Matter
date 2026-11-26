"""Exact U(1)' and Z2 charge selection for chi1 -> chi2 + real phi.

Vectorlike Dirac chi1/chi2, diagonal U(1)' charges and renormalizable
phi * bar(chi2) chi1 + h.c. A single REAL scalar cannot have nonzero
continuous U(1)' charge. This is a scoped analytic model-viability
gate, not a proof against all symmetry-breaking/complex scalar UV models.
"""
from __future__ import annotations
from fractions import Fraction


def analyze(q1,q2,qphi=0,*,real_phi=True):
    q1,q2,qphi=Fraction(q1),Fraction(q2),Fraction(qphi)
    allowed_real=(not real_phi or qphi==0)
    # bar(chi2) transforms with -q2. Charge residual qphi-q2+q1.
    residual=qphi-q2+q1
    invariant=bool(residual==0 and allowed_real)
    chi1_vector=(q1!=0)
    chi2_vector=(q2!=0)
    # Exact Z2: chi1 odd, phi odd, chi2 even. Yukawa parity is +1.
    z2_parity=(-1)*(-1)*(+1)
    return {
      "charges":{"chi1":str(q1),"chi2":str(q2),"phi":str(qphi)},
      "gauge_yukawa_charge_residual":str(residual),
      "one_real_phi_U1_representable":allowed_real,
      "renormalizable_yukawa_gauge_invariant":invariant,
      "z2_yukawa_even":z2_parity==1,
      "chi1_diagonal_vector_current_nonzero":chi1_vector,
      "chi2_diagonal_vector_current_nonzero":chi2_vector,
      "chi1_only_vector_portal_with_neutral_chi2_and_real_phi":
          bool(real_phi and chi1_vector and not chi2_vector and invariant),
      "not_a_full_UV_completion":True
    }


def finite_charge_grid(qmax=3):
    """Finite rational arithmetic control; not the analytic proof."""
    pairs=[]
    for i in range(-qmax,qmax+1):
        for j in range(-qmax,qmax+1):
            a=analyze(i,j)
            if a["renormalizable_yukawa_gauge_invariant"]:
                assert i==j
                if i!=0:
                    assert a["chi2_diagonal_vector_current_nonzero"]
                pairs.append((i,j))
    return pairs


if __name__=="__main__":
    import json
    print(json.dumps({
       "equal_charges":analyze(1,1,0),
       "chi1_only_portal":analyze(1,0,0),
       "complex_charged_phi_escape":analyze(1,0,-1,real_phi=False),
       "grid_invariant_pairs":finite_charge_grid()
    },indent=2))
