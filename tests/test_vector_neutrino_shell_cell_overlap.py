"""Exact angular-cell intersection vs center-hit quadrature, with negative controls."""
import math
import numpy as np
from scipy.special import expit
from kinetics.vector_neutrino_shell_cell_overlap import (
    graph_overlap,shell_angular_normalization,
)
from kinetics.vector_neutrino_reversible_event_graph import (
    synthetic_spectra,event_moments,independent_continuum_weak,
    sweep,ledger,
)


def test_exact_flat_angular_shell_probability_many_resolutions():
    for h in (.4,.2,.1,.05):
        g=graph_overlap(m=.4,h=h,emax=4)
        a=shell_angular_normalization(g)
        assert np.max(abs(a-1))<1e-12
        assert len(g["edges"])>0
    naive=graph_overlap(.4,.2,4,cell_overlap=False)
    assert np.max(abs(shell_angular_normalization(naive)-1))>.1


def test_cell_overlap_exposes_residual_energy_and_cutoff_errors():
    ref=independent_continuum_weak(.4,4,nr=144,nangle=144)
    g=graph_overlap(.4,.2,4)
    n=graph_overlap(.4,.2,4,cell_overlap=False)
    good=event_moments(g,*synthetic_spectra(g))
    bad=event_moments(n,*synthetic_spectra(n))
    assert abs(bad["dNV"]/ref["dNV"]-1)>.05
    assert abs(good["dNV"]/ref["dNV"]-1)<.005
    # Correction is NOT full moment convergence:
    assert abs(good["dEV"]/ref["dEV"]-1)>.01
    fine=graph_overlap(.4,.05,4)
    result=event_moments(fine,*synthetic_spectra(fine))
    assert abs(result["dEV"]/ref["dEV"]-1)<.01
    assert abs(result["dEV"]/ref["dEV"]-1)>.004


def test_midpoint_daughter_graph_is_exact_event_energy_and_positive():
    g=graph_overlap(.4,.2,4)
    E,k=g["E"],g["k"]
    for v,i,l,W in g["edges"]:
        assert math.isclose(k[i]+k[l],E[v],abs_tol=1e-13)
        assert W>0
    b,f,a=synthetic_spectra(g)
    first=ledger(g,b,f,a)
    bn,fn,an=sweep(g,b,f,a,dt=2.)
    last=ledger(g,bn,fn,an)
    for key in ("reduced_charge","lepton_asymmetry","energy"):
        assert math.isclose(first[key],last[key],abs_tol=1e-12,rel_tol=1e-12)
    assert last["entropy"]>first["entropy"]
    assert np.min(bn)>=0 and 0<=np.min(fn)<=np.max(fn)<=1
    assert 0<=np.min(an)<=np.max(an)<=1


def test_chemical_fixed_point_at_nonzero_finite_mu_with_overlap():
    g=graph_overlap(.4,.2,4)
    E,k=g["E"],g["k"]
    f=expit(-k+.1)
    a=expit(-k+.05)
    b=1/np.expm1(E-.15)
    bn,fn,an=sweep(g,b,f,a,dt=1.0)
    assert max(np.max(abs(bn-b)),np.max(abs(fn-f)),
               np.max(abs(an-a)))<1e-12
