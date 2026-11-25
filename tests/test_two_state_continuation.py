"""Checks for new conditional two-state kinetic and mixing bounds."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'kinetics'))
from soft_mixing_gap_tradeoff import threshold_factor,ratios,gcrit
from two_state_halo_kinetic import Box
from mixing_crossing_bound import vacuum_energy

def test_threshold_factor_exact():
    for d in [.001,.1,1.,4.,10.]:
        ref=quad(lambda y:np.sqrt(y*(y+d))*np.exp(-y-d),0.,80.,epsabs=1e-12)[0]
        assert abs(threshold_factor(d)/ref-1)<1e-8
        assert np.exp(-d)<=threshold_factor(d)<=np.exp(-d)*(1+d/2)
    assert abs(gcrit(0.)-2)<1e-12
    assert sum(ratios(.3,3.))<sum(ratios(.3,.1))

def test_two_state_box_local_invariants():
    b=Box(20,.04,.317)
    c=b.validate()
    assert c['balance_relative']<1e-11
    assert c['number_residual']<1e-12
    assert c['energy_residual']<1e-12
    assert c['entropy_production']>0
    assert c['jacobian_relative']<1e-5
    assert c['initial_heavy_gain_per_tau']>0
    e,eq=b.discrete_eq()
    assert eq is not None and .3<eq['heavy_fraction']<.6
    assert b.vol@abs(b.f(e,compiled=False))<1e-12

def test_rationalized_vacuum_gap():
    D=2.5e-6;d=3.2e-6
    assert np.isclose(vacuum_energy(0.,0.,d,D),D/(2*86),rtol=1e-8)
    assert np.isclose(vacuum_energy(1e8,1.,d,D),np.hypot(D,2*d)/(2e8),rtol=1e-11)
