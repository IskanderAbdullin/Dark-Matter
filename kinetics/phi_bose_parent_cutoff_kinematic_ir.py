"""Sharp finite-parent-cutoff infrared support for boosted on-shell 1->2 Bose daughters.

This is a strictly kinematic source-support result, not BEC or a physical
freezeout. The positive massive daughter has a finite energy floor when
parent momentum is restricted; physical Bose occupation exists below that
floor but the truncated decay operator cannot reach it.
"""
import json
import math

M1, M2, MPH = .1, .02, 86e-9


def daughter_energy_bounds(p, m1=M1, m2=M2, mphi=MPH):
    if not (p >= 0 and m1 > m2 + mphi and m2 >= 0 and mphi > 0):
        raise ValueError('physical open two-body decay with positive boson mass required')
    kallen = (m1*m1 - (m2+mphi)**2)*(m1*m1 - (m2-mphi)**2)
    qstar = math.sqrt(kallen)/(2*m1)
    estar = math.hypot(qstar, mphi)
    eparent = math.hypot(p, m1)
    a = eparent*estar/m1
    b = p*qstar/m1
    # Rationalization avoids loss of precision at large p or tiny mphi.
    emin = (p*p*mphi*mphi + m1*m1*estar*estar)/(m1*(eparent*estar + p*qstar))
    emax = a+b
    turning_momentum = m1*qstar/mphi
    return dict(p=p, minimum=emin, maximum=emax, qstar=qstar,
                Ephi_star=estar, turning_momentum=turning_momentum)


def min_energy_with_parent_cutoff(pmax, m1=M1, m2=M2, mphi=MPH):
    if pmax <= 0: raise ValueError('positive momentum cutoff required')
    d=daughter_energy_bounds(pmax,m1,m2,mphi)
    # dEmin/dp = (p/E1 Ephi* - p*)/m1. Sign changes exactly at
    # p_turn = m1 p*/mphi, where Emin reaches the boson rest mass.
    if pmax >= d['turning_momentum']:
        return mphi
    return d['minimum']


def quadratic_first_cell(n, pmax, m1=M1, m2=M2, mphi=MPH):
    """Use the identical scalar mesh definition from phi_bose_cell_entropy.fixture."""
    if n < 2 or pmax <= 0:raise ValueError('nontrivial finite scalar grid')
    upper=math.hypot(pmax,m1)+.2
    edge=mphi+(upper-mphi)/n**2
    support_min=min_energy_with_parent_cutoff(pmax,m1,m2,mphi)
    critical_n=math.sqrt((upper-mphi)/(support_min-mphi)) if support_min>mphi else math.inf
    return {'n':n,'pmax_GeV':pmax,'first_cell_upper_GeV':edge,
            'true_min_daughter_energy_GeV':support_min,
            'first_cell_has_no_possible_on_shell_event':bool(edge < support_min),
            'critical_n':critical_n}


def audit():
    return {'schema_version':1,'model':{'m1_GeV':M1,'m2_GeV':M2,'mphi_GeV':MPH},
            'on_shell_extrema':[daughter_energy_bounds(p) for p in (0.,.1,1.,2.,8.,16.)],
            'parent_cutoffs':[{'pmax_GeV':p,'kinematic_minimum_GeV':min_energy_with_parent_cutoff(p),
                              'mesh_cells':[quadratic_first_cell(n,p) for n in (16,24,32,40,42,43,48,64,128,256,512)]}
                             for p in (2.,8.,16.)],
            'interpretation':'At pmax2, N48 and N64 first cells are kinematically inaccessible to any on-shell source: exactly zero support independently of angular quadrature, not an empty-sampled cell or evidence of BEC. For pmax8/16 the thresholds move.',
            'class_ready':False}


if __name__=='__main__':
    print(json.dumps(audit(),indent=2,allow_nan=False))
