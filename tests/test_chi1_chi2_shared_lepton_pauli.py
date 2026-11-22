"""Strong paired χ1+χ2 time is nonterminal; both use one lepton f."""
import math
import numpy as np
from scipy.special import expit
from kinetics.chi1_chi2_shared_lepton_pauli import (
    system,ledger,rhs,evolve,interp_fermi,cell_weights
)
from kinetics.independent_chi_lepton_5d_quadrature import evaluate


def test_shared_lepton_two_species_strong_weak_heat_and_entropy():
    cfg=system(bins=10)
    chi=cfg['initial_chi'];lep=cfg['initial_lepton']
    dc,dl,r=rhs(cfg,chi,lep)
    assert len(dc)==2
    assert all(abs(x)<1e-24 for x in r['Nchi_rates'])
    assert abs(r['Nlep_rate'])<1e-24
    assert abs(sum(r['Qchi'])+r['Qlep'])<1e-24
    assert r['entropy_direct']>0
    assert abs(r['entropy_node']-r['entropy_direct'])<2e-25
    h=interp_fermi(cfg['lg'],lep)
    separate=[]
    for M,grid,fi in zip(cfg['masses'],cfg['chi_grid'],chi):
        f=interp_fermi(grid,fi)
        separate.append(evaluate((f,h),M=M,m=cfg['mlepton'],nr=4,nmu=4,
             nc=8,nphi=8,cutoff=12.,angular_scheme='log_transfer'))
    for i in (0,1):
        assert math.isclose(separate[i]['Qchi'],r['Qchi'][i],rel_tol=2e-11)
    assert math.isclose(sum(x['entropy'] for x in separate),
                        r['entropy_direct'],rel_tol=2e-11)


def test_common_temperature_with_three_independent_mu_is_fixed_point():
    cfg=system(bins=10,mlepton=.5)
    chi=[expit(.1-g) for g in cfg['chi_grid'][:1]]+[
         expit(-.05-cfg['chi_grid'][1])]
    lep=expit(-.08-cfg['lg'])
    dc,dl,r=rhs(cfg,chi,lep)
    assert max(max(abs(z)) for z in dc)<1e-20
    assert max(abs(dl))<1e-20
    assert abs(r['entropy_direct'])<1e-32


def test_pauli_bounded_nonlinear_time_all_three_numbers_energy_entropy():
    cfg=system(bins=10)
    out=evolve(cfg,total_time=1e10,atol=1e-7)
    assert out['accepted']>=5 and out['rejected']>=1
    assert 0<out['min_step_entropy_gain'] and 0<out['entropy_gain']
    assert out['max_local_spectral_error']<=1e-7
    for key in ('Nchi1_drift','Nchi2_drift','Nlepton_drift','energy_drift'):
        assert abs(out[key])<2e-12
    assert min(out['final_fchi1_min'],out['final_fchi2_min'],
               out['final_lepton_min'])>0


def test_physical_separate_radial_capacity_normalized_to_density_volume():
    cfg=system(bins=10,mlepton=.5)
    for M,w,gr in zip(cfg['masses'],cfg['chi_w'],cfg['chi_grid']):
        expected=4/(6*math.pi**2)*(gr[-1]**2-M*M)**1.5
        assert math.isclose(sum(w),expected,rel_tol=2e-14)
    m=cfg['mlepton']
    expected=4/(6*math.pi**2)*(cfg['lg'][-1]**2-m*m)**1.5
    assert math.isclose(sum(cfg['lw']),expected,rel_tol=2e-14)
