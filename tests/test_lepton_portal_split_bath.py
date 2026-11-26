"""Two-bath lepton vector-portal Pauli heat: sign-flip negative control."""
from kinetics.lepton_portal_split_bath import (
  TNU_OVER_TGAMMA_INSTANTANEOUS,split_bath_heat)

def test_nonisothermal_neutrino_bath_can_reverse_net_dark_heat_sign():
    a=split_bath_heat(.1,n=20)
    assert abs(a["Tnu_over_Tgamma"]-TNU_OVER_TGAMMA_INSTANTANEOUS)<1e-14
    assert a["Q_e_over_G2_me9"]>0
    assert a["Q_two_neutrino_flavors_over_G2_me9"]<0
    assert a["Q_total_over_G2_me9"]<0
    assert a["Qnu_over_Qe"]<-500

def test_same_bath_temperature_restores_positive_neutrino_heat():
    a=split_bath_heat(.1,Tnu_over_Tgamma=1.,n=20)
    assert a["Q_e_over_G2_me9"]>0
    assert a["Q_two_neutrino_flavors_over_G2_me9"]>0
