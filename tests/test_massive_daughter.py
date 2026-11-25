"""Deterministic finite-mass and conditional PSD checks; no viability claim."""
import json
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'kinetics'))
sys.path.insert(0,str(ROOT/'cosmology'))

from two_mass_bath import rectangle_areas
from two_mass_geometry import cell_area
from yukawa_rate_bath import YukawaRateBath
from thermal_consistency_screen import backward_energy
from daughter_bath_class import canonical,free_streaming
from three_scale_psd import export_three
from kinetic_psd_class import transformed_quadrature_check


def test_near_degenerate_geometry():
    x,m,b=1.,.01,.98
    e=np.geomspace(1e-5,8,40);r=np.geomspace(2e-5,9,42)
    a=rectangle_areas(e,r,x,m,b)
    rev=rectangle_areas(r,e,x,b,m).T
    scale=np.diff(np.hypot(e,m))[:,None]*np.diff(np.hypot(r,b))[None,:]
    assert np.min(a)>=0
    assert np.max(np.abs(a-rev)/scale)<5e-5
    idx=np.argwhere((a>0)&(a<.99*scale))
    assert len(idx)>0
    for i,j in idx[::max(1,len(idx)//15)]:
        ref=cell_area(np.hypot(e[i],m),np.hypot(e[i+1],m),np.hypot(r[j],b),np.hypot(r[j+1],b),x,m,b)
        assert abs(a[i,j]-ref)/scale[i,j]<5e-5


def test_empty_parent_early_source():
    model=YukawaRateBath(n=40,scalar_ratio=1e-12,fermion_ratio=.9,kappa=.00001,
        kappa2=.03,kappa_sm1=47.,gamma=.01,scatter=True,qmin=1e-7)
    r=model.run(.1,.15,2e-6,1e-15,False,86000.,'empty')
    assert r['initial_mother_number']==0
    f=r['final']
    assert f['mother_number']>0 and f['number']>0
    assert abs(f['number_source_residual'])<1e-6
    assert abs(f['continuity_residual'])<5e-5


def test_backward_energy_stable_formula():
    assert np.isclose(backward_energy(1.,.9,0.,3.),.095/(np.hypot(1,3)+3),rtol=1e-14)
    assert backward_energy(1.,.9,.01,3.)>=.01


def test_three_scale_export_actual_candidate(tmp_path):
    rows=json.loads((ROOT/'kinetics/results/yukawa-degeneracy-2026-11-25.json').read_text())['outputs']['yukawa90_320.json']
    row=canonical(rows[0])
    e=export_three(row,tmp_path)
    err=transformed_quadrature_check(e,120)
    assert len(e['components'])==3
    assert np.max(np.abs(err))<.001
    assert max(map(abs,e['moment_relative_errors']))<.006
    fs=free_streaming(row)
    assert .001<fs['through_equality']['quantiles_Mpc']['0.5']<.003
