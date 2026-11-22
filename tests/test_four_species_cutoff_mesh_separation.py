"""Separates moving-grid artefacts from physical cutoff for synthetic χ1χ2φℓ."""
import math
import pytest
from kinetics.four_species_rate_error_budget import sample
from kinetics.four_species_cutoff_mesh_separation import compare_mesh,independent_weak


def test_affine_fd_logit_is_mesh_invariant_at_fixed_incoming_phase_space():
    a=sample(cut=8.,bins=32,nlocal=3,nmu=6,nc=12,nphi=12,shape='FD')
    b=sample(cut=8.,mesh_cut=20.,bins=32,nlocal=3,nmu=6,nc=12,nphi=12,shape='FD')
    for field in ('Qchi1_elastic','Qchi2_elastic','entropy_elastic','Yukawa_Rphi'):
        assert math.isclose(a[field],b[field],rel_tol=2e-11)


def test_nonthermal_rate_and_entropy_are_confounded_by_moving_mesh():
    a=sample(cut=12.,bins=64,nlocal=3,nmu=6,nc=12,nphi=12,shape='oscillatory')
    b=sample(cut=12.,mesh_cut=20.,bins=64,nlocal=3,nmu=6,nc=12,nphi=12,shape='oscillatory')
    assert abs(a['Qchi1_elastic']/b['Qchi1_elastic']-1)>.05
    assert abs(a['entropy_elastic']/b['entropy_elastic']-1)>.2
    assert max(abs(a['combined_energy_residual']),abs(b['combined_energy_residual']))<5e-24
    assert a['entropy_elastic']>0 and b['entropy_elastic']>0


def test_four_grid_entropy_is_not_converged_while_5d_weak_is_independent():
    x=compare_mesh(bins=(32,64,128,256))
    rows=x['projected_5D_mode_occupations']
    errors=[abs(r['entropy_error_relative']) for r in rows]
    assert all(a>b for a,b in zip(errors,errors[1:]))
    assert errors[0]>.8
    assert errors[-1]>.04  # even 256-mode entropy still biased >4%!
    assert abs(rows[-1]['Qchi1_error_relative'])<.01
    assert abs(rows[-1]['Qchi2_error_relative'])<.01
    assert all(r['entropy_identity_error']<5e-23 for r in rows)
    assert all(abs(r['energy_residual'])<5e-23 for r in rows)


def test_fail_closed_if_collision_cutoff_exceeds_shared_energy_mesh():
    with pytest.raises(ValueError):
        sample(cut=12.,mesh_cut=10.,shape='FD',nlocal=3,nmu=6,nc=12,nphi=12)
