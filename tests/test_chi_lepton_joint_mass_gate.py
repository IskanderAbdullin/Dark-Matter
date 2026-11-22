"""Both χ masses and finite lepton mass are conditional, not Γelastic."""
from kinetics.chi_lepton_joint_mass_gate import case
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def test_two_chi_species_and_lepton_mass_conserve_pair_energy_and_entropy():
    for M in (.1,.02):
        for m in (0.,.5):
            x=case(M,m,"smooth","coarse")
            assert x["Qchi"]>0 and x["Qlepton"]==-x["Qchi"]
            assert x["entropy"]>0 and x["negative_entropy_events"]==0
            assert x["max_energy_defect"]<1e-13


def test_mass_effect_is_separate_from_chemical_detailed_balance():
    for M in (.1,.02):
        eq=evaluate("fd_eq_chem",M=M,m=.5,mv=.1,
                    nr=4,nmu=6,nc=12,nphi=12,
                    radial_scheme="segmented_p",segments=8,
                    angular_scheme="log_transfer")
        assert abs(eq["Qchi"])<1e-26 and abs(eq["entropy"])<1e-34
    massless=case(.1,0.,"smooth","coarse")["Qchi"]
    massive=case(.1,.5,"smooth","coarse")["Qchi"]
    assert .85<massive/massless<.95


def test_multiaxis_changes_do_not_prove_physical_continuum_error():
    a=case(.02,.5,"oscillatory","coarse")
    b=case(.02,.5,"oscillatory","mid")
    assert a["Qchi"]>0 and b["Qchi"]>0
    assert abs(b["Qchi"]/a["Qchi"]-1)>.02
