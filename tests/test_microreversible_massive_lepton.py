"""Independent zero-mass, reciprocity, threshold and finite-T Pauli controls."""
import math,sys
from pathlib import Path
import numpy as np
import pytest
from scipy.special import expit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from microreversible_massive_lepton import (
    massive_vector_kernel,massive_chie_lattice,thermal_heat_example,ELECTRON_MASS_GEV)
from microreversible_chie_lattice import invariant_kernel,sources


def test_massless_limit_and_direct_time_reversal_massive():
    for E,k,Ep,kp in [(1.,2.,2.,1.),(.7,1.2,1.4,.5),(3.,4.,4.,3.)]:
        a=massive_vector_kernel(.1,0.,.001,E,k,np.array([Ep]),np.array([kp]))
        b=invariant_kernel(.1,.001,E,k,np.array([Ep]),np.array([kp]))
        np.testing.assert_allclose(a,b,rtol=2e-13,atol=0)
        m=ELECTRON_MASS_GEV
        fw=massive_vector_kernel(.1,m,.001,E,k,np.array([Ep]),np.array([kp]))[0]
        rev=massive_vector_kernel(.1,m,.001,Ep,kp,np.array([E]),np.array([k]))[0]
        assert fw>0 and math.isclose(fw,rev,rel_tol=1e-12)


def test_common_T_fd_stops_for_independent_mu():
    g=massive_chie_lattice(T=.001,dE=.00025,nchi=50,ne=50)
    E=g["Echi"];k=g["Ee"];T=.001
    for muc,mue in [(0.,0.),(.009,-.001)]:
        r=sources(g,expit((muc-E)/T),expit((mue-k)/T))
        assert abs(r["Qchi_GeV5"])<1e-40
        assert abs(r["Nchi_source_GeV4"])<1e-40
        assert abs(r["Ne_source_GeV4"])<1e-40
        assert np.all(g["Ve"]>0)


def test_electron_mass_suppresses_charge_symmetric_thermal_exchange():
    T=.0001
    m0=thermal_heat_example(T=T,m=0.,n=80)
    m1=thermal_heat_example(T=T,m=ELECTRON_MASS_GEV,n=80)
    assert .2<m1["Qchi"]/m0["Qchi"]<.5
    assert m0["Qchi"]>m1["Qchi"]>0
    assert m1["Qe"]==-m1["Qchi"]
    assert m1["entropy"]>0


def test_low_T_suppression_stronger_than_at_T_greater_than_mass():
    hi=thermal_heat_example(T=.0002,n=60)
    hi0=thermal_heat_example(T=.0002,m=0.,n=60)
    lo=thermal_heat_example(T=.00005,n=60)
    lo0=thermal_heat_example(T=.00005,m=0.,n=60)
    assert .6<hi["Qchi"]/hi0["Qchi"]<1.
    assert 0<lo["Qchi"]/lo0["Qchi"]<.06
    assert lo["Qchi"]/lo0["Qchi"]<hi["Qchi"]/hi0["Qchi"]


@pytest.mark.parametrize("kw",[{"M":0.},{"m":-1.},{"G":-1.},{"nchi":2},
                                {"ne":2},{"dE":0.}])
def test_invalid_grid_rejected(kw):
    with pytest.raises(ValueError):massive_chie_lattice(**kw)
