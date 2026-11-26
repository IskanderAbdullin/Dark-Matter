"""Exact boosted 1->2 daughter infrared accessibility threshold.

This is kinematics only. It does not calculate an evolving distribution,
a Bose condensate, a scattering cross section, or a relic abundance.
All energies/momenta use GeV.
"""
import math

M1 = 0.1
M2 = 0.02
MPHI = 86e-9


def rest_kinematics(m1=M1, m2=M2, mphi=MPHI):
    if not (m1 > m2 + mphi > 0 and mphi > 0):
        raise ValueError("Physical open decay with massive daughter required")
    estar = (m1*m1 + mphi*mphi - m2*m2)/(2*m1)
    qstar = math.sqrt((estar-mphi)*(estar+mphi))
    return estar, qstar


def parent_threshold_for_daughter_energy(k, m1=M1, m2=M2, mphi=MPHI):
    """Smallest parent p allowing ANY scalar lab energy <= k.

    For k>=E_phi*, parent at rest already produces an allowed scalar.
    For k<mphi no on-shell scalar exists. Includes the exact scalar-rest
    threshold at k=mphi (extremely large parent momentum for tiny mphi).
    """
    estar, qstar = rest_kinematics(m1,m2,mphi)
    if k < mphi:
        raise ValueError("Below massive scalar on-shell rest energy")
    if k >= estar:
        return 0.0
    return m1*math.sinh(math.asinh(qstar/mphi)-math.acosh(k/mphi))


def daughter_min_at_parent(p, m1=M1, m2=M2, mphi=MPHI):
    """Minimum scalar lab energy at a fixed parent momentum.

    Numerically stable positive rationalization prevents catastrophic
    cancellation of E*E1 - p*q* at ultra-relativistic parent boosts.
    """
    if p < 0:
        raise ValueError("Parent momentum must be nonnegative")
    estar, qstar = rest_kinematics(m1,m2,mphi)
    e1 = math.hypot(p,m1)
    return (mphi*mphi*p*p + m1*m1*estar*estar)/(m1*(e1*estar+p*qstar))


def min_over_parent_cutoff(pmax, m1=M1, m2=M2, mphi=MPHI):
    """Minimum energy over 0<=p<=pmax, including the nonmonotone branch."""
    if pmax < 0:
        raise ValueError("Cutoff must be nonnegative")
    pstar=parent_threshold_for_daughter_energy(mphi,m1,m2,mphi)
    return mphi if pmax>=pstar else daughter_min_at_parent(pmax,m1,m2,mphi)


def table():
    thresholds=[MPHI,2*MPHI,10*MPHI,1e-6,1e-5,1e-4,1e-3]
    return {"thresholds": [{"scalar_energy_GeV":k,
                             "minimum_parent_p_GeV":parent_threshold_for_daughter_energy(k)}
                            for k in thresholds],
            "cutoffs": [{"parent_pmax_GeV":p,
                         "minimum_scalar_energy_GeV":min_over_parent_cutoff(p)}
                        for p in (2.,8.,16.,24.)]}


if __name__ == "__main__":
    import json
    print(json.dumps(table(),indent=2))
