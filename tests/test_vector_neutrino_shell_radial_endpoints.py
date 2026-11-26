"""Vector finite-cutoff radial quadrature: genuine endpoint bias, not fixed Γel."""
import numpy as np
from kinetics.vector_neutrino_shell_radial_endpoints import adjusted_graph,scan
from kinetics.vector_neutrino_shell_cell_overlap import shell_angular_normalization
from kinetics.vector_neutrino_reversible_event_graph import synthetic_spectra,sweep,ledger


def test_radial_endpoint_not_same_as_full_bin_and_ang_norm_stays_exact():
    for scheme in ("original","upper_trapezoid","shell_cell_integral"):
        g=adjusted_graph(.4,.2,4,scheme=scheme)
        assert np.max(abs(shell_angular_normalization(g)-1))<2e-12
        assert np.min(g["wV"])>0
        assert all(W>0 for _,_,_,W in g["edges"])
    orig=adjusted_graph(.4,.2,4,scheme="original")
    half=adjusted_graph(.4,.2,4,scheme="upper_trapezoid")
    assert half["wV"][-1]==.5*orig["wV"][-1]


def test_endpoint_correction_reduces_energy_error_not_all_moments():
    rows=scan(.4,4,steps=(.2,.1,.05))["scan"]
    for row in rows:
        orig=abs(row["original"]["dEV_rel"])
        trap=abs(row["upper_trapezoid"]["dEV_rel"])
        exact=abs(row["shell_cell_integral"]["dEV_rel"])
        assert exact<trap<orig
    fine=rows[-1]
    assert abs(fine["shell_cell_integral"]["dEV_rel"])<4e-4
    # Even the better cell-weighting is NOT continuum-exact in number:
    assert abs(fine["shell_cell_integral"]["dNV_rel"])>1e-3


def test_exact_event_invariants_and_entropy_survive_radial_change():
    for scheme in ("upper_trapezoid","shell_cell_integral"):
        g=adjusted_graph(.4,.2,4,scheme=scheme)
        b,f,a=synthetic_spectra(g)
        before=ledger(g,b,f,a)
        B,F,A=sweep(g,b,f,a,dt=2.)
        after=ledger(g,B,F,A)
        assert after["entropy"]>before["entropy"]
        for metric in ("energy","reduced_charge","lepton_asymmetry"):
            assert abs(after[metric]-before[metric])<1e-12
        assert np.min(B)>=0 and 0<=np.min(F)<=np.max(F)<=1
        assert 0<=np.min(A)<=np.max(A)<=1
