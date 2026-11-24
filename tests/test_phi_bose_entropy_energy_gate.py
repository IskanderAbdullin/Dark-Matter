"""Targeted regression for conditional ideal-gas entropy bound."""
from kinetics.phi_bose_entropy_energy_gate import main

def test_entropy_energy_gate():
    r = main()
    assert all(x['trial_bound'] <= x['rho_over_one_nu']+1e-8 for x in r['bbn'])
    assert abs(r['g_dark_dec']-7.991670794968991)<3e-9
    assert abs(r['bbn'][-2]['rho_over_one_nu']-.88751172857)<1e-8
    assert r['today_bound']['omega_h2_min']<.12<r['relativistic_scalar_only_freezeout']['omega_h2']
    assert r['count_constrained_gate']['rho_lower_over_one_nu']>47
    assert abs(r['bbn'][-2]['rho_over_one_nu']/(4/7*(10.75*r['s_dark_over_s_sm'])**(4/3))-1)<1e-6
