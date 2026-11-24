"""Independent exact massive daughter inverse parent rapidity support tests."""
import math
import numpy as np
from kinetics.phi_bose_invariant_parent_support import (
    parent_interval,global_support_floor,
)
from kinetics.phi_bose_physical_cutcell_weak import (
    M1,M2,MPHI,KSTAR,shell,
)


def _boosted_energy_interval(p,mass):
    e1=math.hypot(p,M1)
    estar=math.hypot(KSTAR,mass)
    return ((e1*estar-p*KSTAR)/M1,
            (e1*estar+p*KSTAR)/M1)


def test_general_inverse_rapidity_support_both_daughters():
    for mass in (MPHI,M2):
        for P in (2.,8.,16.):
            true_floor=global_support_floor(mass,P)
            pstar=M1*KSTAR/mass
            pbest=min(P,pstar)
            assert math.isclose(true_floor,
                                _boosted_energy_interval(pbest,mass)[0],
                                rel_tol=5e-9)
            for E in (max(mass*(1+1e-6),true_floor*1.01),
                      .01,.05,.1,.2,1.):
                if E<=mass:continue
                bounds=parent_interval(E,mass,P)
                for p in np.linspace(0,P,25):
                    lo,hi=_boosted_energy_interval(p,mass)
                    valid=lo-1e-11<=E<=hi+1e-11
                    by_bounds=(bounds is not None and
                               bounds[0]-1e-11<=p<=bounds[1]+1e-11)
                    assert valid==by_bounds,(mass,P,E,p,bounds,lo,hi)


def test_scalar_rest_requires_55p8TeV_and_upper_parent_bound_beyond():
    pstar=M1*KSTAR/MPHI
    assert 55800<pstar<55830
    for P in (2.,8.,16.):
        assert parent_interval(MPHI,MPHI,P) is None
        assert math.isclose(global_support_floor(MPHI,P),
                            shell(P)[1],rel_tol=1e-10)
    bounds=parent_interval(86.1e-9,MPHI,60000.)
    assert 53000<bounds[0]<54000
    assert 58000<bounds[1]<59000
    assert bounds[1]<60000
    for P in (60000.,100000.):
        assert math.isclose(global_support_floor(MPHI,P),
                            MPHI,rel_tol=1e-12)
        assert shell(P)[1]>global_support_floor(MPHI,P)


def test_fermionic_daughter_rest_at_0p24GeV_has_two_parent_limits():
    pstar=M1*KSTAR/M2
    assert math.isclose(pstar,.24,rel_tol=1e-11)
    bounds=parent_interval(.0201,M2,2.)
    assert .21<bounds[0]<.22 and .26<bounds[1]<.27
    E2star=(M1*M1+M2*M2-MPHI*MPHI)/(2*M1)
    other=parent_interval(E2star,M2,2.)
    assert other[0]<1e-10 and 1.24<other[1]<1.26
