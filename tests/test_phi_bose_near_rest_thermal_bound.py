"""Fixed thermal FD/BE near-rest Yukawa event bound regression tests."""
import math
from kinetics.phi_bose_near_rest_thermal_bound import thermal_IR_upper_bound
from kinetics.phi_bose_invariant_parent_support import parent_interval
from kinetics.phi_bose_physical_cutcell_weak import MPHI,M1,KSTAR


def test_nested_exact_parent_intervals_and_physical_boost():
    for top in (MPHI+1e-10,MPHI+1e-9,MPHI+1e-8,MPHI+1e-7):
        lower,upper=parent_interval(top,MPHI,1e10)
        assert 0<lower<upper<1e7
        for k in (MPHI+1e-14,(top+MPHI)/2,top):
            lo,hi=parent_interval(k,MPHI,1e10)
            assert lower<=lo+1e-5 and hi<=upper+1e-5
            p=(lo+hi)/2
            e1=math.hypot(p,M1)
            estar=math.hypot(KSTAR,MPHI)
            A=e1*estar/M1;B=p*KSTAR/M1
            assert A-B-1e-7<=k<=A+B+1e-7


def test_log_bound_enormously_below_relevant_thermal_weak_rates():
    row=thermal_IR_upper_bound(MPHI+1e-10)
    assert 53186<row['parent_allowed_union_GeV'][0]<53187
    assert 58571<row['parent_allowed_union_GeV'][1]<58572
    assert math.isclose(row['log10_absolute_rate_bound_GeV4'],
                        -23093.306185861147,abs_tol=1e-7)
    assert row['log10_absolute_scalar_heat_bound_GeV5']<row['log10_absolute_rate_bound_GeV4']-7
    for top in (MPHI+1e-9,MPHI+1e-8,MPHI+1e-7):
        later=thermal_IR_upper_bound(top)
        assert row['log10_absolute_rate_bound_GeV4']<later['log10_absolute_rate_bound_GeV4']
        row=later


def test_nearly_critical_still_finite_bose_upper_bound():
    regular=thermal_IR_upper_bound(MPHI+1e-10)
    nearly=thermal_IR_upper_bound(MPHI+1e-10,muphi=MPHI-1e-16)
    assert nearly['bose_thermal_max']>1e15
    assert nearly['log10_absolute_rate_bound_GeV4']>regular['log10_absolute_rate_bound_GeV4']
    assert nearly['log10_absolute_rate_bound_GeV4']<-23000


def test_bound_explicitly_excludes_critical_be_chemical_potential():
    try: thermal_IR_upper_bound(MPHI+1e-10,muphi=MPHI)
    except AssertionError:pass
    else:raise AssertionError('A critical φ chemical potential is outside this bound')
