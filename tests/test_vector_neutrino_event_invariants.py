"""Exact modular stoichiometric invariants: flag frozen ν endpoint bins."""
import math
import numpy as np
import pytest
from kinetics.vector_neutrino_event_invariants import (
 audit,stoichiometry,exact_mod_prime_rank,
)
from kinetics.vector_neutrino_reversible_event_graph import graph
from kinetics.vector_neutrino_shell_cell_overlap import graph_overlap


@pytest.mark.parametrize("h",(.4,.2,.1))
def test_integer_endpoint_has_two_spurious_exact_conserved_modes(h):
    x=audit(m=.4,h=h,emax=4,shifted=False)
    y=audit(m=.4,h=h,emax=4,shifted=True)
    assert x["nullity"]==5
    assert y["nullity"]==3
    assert x["extra_exact_invariants_beyond_C_L_E"]==2
    assert y["extra_exact_invariants_beyond_C_L_E"]==0
    assert x["isolated_nu_indices"]==[x["nu_modes_each"]-1]
    assert x["isolated_antinu_indices"]==[x["nu_modes_each"]-1]
    assert y["isolated_nu_indices"]==[]
    assert y["isolated_antinu_indices"]==[]
    assert all(abs(v)<1e-13 for v in x["physical_invariant_residuals"].values())
    assert all(abs(v)<1e-13 for v in y["physical_invariant_residuals"].values())
    for z in (x,y):
        assert z["rank_mod_1000003"]==z["rank_mod_1000033"]
        assert z["rank_mod_1000003"]==z["rank_float_svd"]


def test_stoichiometric_zero_endpoint_row_is_an_independent_invariant():
    g=graph(.4,.2,4)
    S=stoichiometry(g)
    nv=len(g["E"]);nk=len(g["k"])
    assert np.all(S[nv+nk-1]==0)
    assert np.all(S[nv+2*nk-1]==0)
    shifted=stoichiometry(graph_overlap(.4,.2,4))
    assert not np.any(np.all(shifted==0,axis=1))


def test_rank_modular_is_not_floating_point_threshold_guess():
    S=np.array([[-1,0,1],[1,1,0],[0,-1,-1]],dtype=np.int64)
    for prime in (1000003,1000033):
        assert exact_mod_prime_rank(S,prime)==2
