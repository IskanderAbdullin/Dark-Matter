"""A narrow IR Pauli packet requires adapted quadrature; upper bound not eigenmode."""
from __future__ import annotations
import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"kinetics"))
from localized_pauli_packets import packet_rayleigh


def test_split_IR_quadrature_refinement_and_tail_cutoff():
    a=packet_rayleigh(.25,.08,n=20,cutoff=32,split=True)
    b=packet_rayleigh(.25,.08,n=24,cutoff=32,split=True)
    c=packet_rayleigh(.25,.08,n=20,cutoff=40,split=True)
    assert a["Rayleigh_GeV"]>0 and b["Rayleigh_GeV"]>0
    assert abs(a["Rayleigh_GeV"]/b["Rayleigh_GeV"]-1)<.002
    assert abs(c["Rayleigh_GeV"]/a["Rayleigh_GeV"]-1)<.001
    assert b["not_full_spectral_gap"]


def test_IR_packet_upper_bound_below_prior_11_mode_Ritz():
    slow=packet_rayleigh(.25,.08,n=20,cutoff=32,split=True)
    broad=packet_rayleigh(.5,.5,n=20,cutoff=32,split=True)
    assert slow["Rayleigh_GeV"]<.65*8.564757722572995e-8
    assert slow["Rayleigh_GeV"]<broad["Rayleigh_GeV"]


def test_bad_uniform_initial_p_quadrature_fakes_too_small_kinetic_rate():
    coarse=packet_rayleigh(.25,.08,n=20,cutoff=32,split=False)
    good=packet_rayleigh(.25,.08,n=20,cutoff=32,split=True)
    assert coarse["Rayleigh_GeV"]<.7*good["Rayleigh_GeV"]
    assert coarse["S_GeV3"]==good["S_GeV3"]


def test_coupling_scaling_and_invalid_width():
    a=packet_rayleigh(.25,.08,G=.001,n=12,split=True)
    b=packet_rayleigh(.25,.08,G=.002,n=12,split=True)
    assert math.isclose(b["Rayleigh_GeV"]/a["Rayleigh_GeV"],4.,rel_tol=2e-12)
    with pytest.raises(ValueError):packet_rayleigh(.2,0)
    with pytest.raises(ValueError):packet_rayleigh(.2,-1)
