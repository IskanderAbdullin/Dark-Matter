"""Quadrature convergence and exact physical radial volume are separate gates."""
import math
from kinetics.chi_lepton_segmented_5d_gate import compare,radial_volume_error
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def test_segmented_p_gauss_preserves_radial_number_measure():
    for m in (0.,.1,.5):
        assert abs(radial_volume_error(m,12.,12,5,"segmented_p"))<3e-14
    # Energy Gauss sees the non-polynomial m>0 threshold square root.
    assert abs(radial_volume_error(.5,12.,12,5,"segmented_energy"))>1e-7


def test_oscillatory_segments_remove_wrong_underresolved_uniform_sign():
    under=evaluate("oscillatory",nr=24,nmu=8,nc=32,nphi=16)
    segmented=compare("oscillatory",nlocal=5)
    assert under["Qchi"]<0<segmented["Qchi"]
    assert segmented["radial_points_per_incoming_species"]==60
    assert math.isclose(segmented["Qchi"],4.568755689044902e-15,rel_tol=3e-10)
    assert segmented["entropy"]>0 and segmented["negative_entropy_events"]==0
    assert segmented["max_energy_defect"]<1e-13


def test_fixed_energy_segment_width_in_cutoff_comparison():
    # Cutoff controls must not secretly change segment energy resolution.
    a=compare("oscillatory",cutoff=12.)
    b=compare("oscillatory",cutoff=16.)
    c=compare("oscillatory",cutoff=20.)
    assert [x["radial_points_per_incoming_species"] for x in (a,b,c)]==[60,80,100]
    assert all(x["Qchi"]>0 and x["entropy"]>0 for x in (a,b,c))
    assert abs(c["Qchi"]/b["Qchi"]-1)<.003
    assert abs(c["entropy"]/b["entropy"]-1)<.003
    # This is NOT a coupled outgoing-angle convergence certificate.


def test_independent_chemical_balance_with_segmentation():
    eq=evaluate("fd_eq_chem",nr=5,nmu=8,nc=32,nphi=16,
                radial_scheme="segmented_p",segments=8)
    assert abs(eq["Qchi"])<1e-27
    assert abs(eq["entropy"])<1e-34
