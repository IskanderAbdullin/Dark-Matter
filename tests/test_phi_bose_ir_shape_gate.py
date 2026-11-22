"""Independent synthetic physical Yukawa φ IR spectrum and UV-cutoff gates."""
import math
from kinetics.phi_bose_ir_shape_gate import rate,reference,cutoff_audit


def test_common_fd_affine_inverse_fugacity_is_exact_on_16_64_256_bins():
    ref=reference('thermal',nrad=24,nangle=48)
    for n in (16,64,256):
        r=rate('thermal',bins=n,nrad=24,nangle=48)
        assert math.isclose(r['Rphi'],ref['R_phi'],rel_tol=1e-12)
        assert math.isclose(r['H'],ref['entropy'],rel_tol=1e-12)
        assert abs(r['Q_residual'])<1e-24


def test_legal_curved_phi_inverse_fugacity_is_grossly_underresolved_on_coarse_bins():
    ref=reference('infrared_curved',nrad=24,nangle=48)
    q=[rate('infrared_curved',bins=n,nrad=24,nangle=48)
       for n in (16,32,64,128,256,512)]
    assert ref['R_phi']<0 and ref['entropy']>0
    assert all(x['H']>0 and x['Rphi']<0 for x in q)
    assert abs(q[0]['Rphi']/ref['R_phi']-1)>.6
    assert abs(q[1]['Rphi']/ref['R_phi']-1)>.25
    assert abs(q[-1]['Rphi']/ref['R_phi']-1)<.01
    assert abs(q[-1]['H']/ref['entropy']-1)<.01
    assert abs(q[-1]['Q_residual'])<1e-23


def test_broad_phi_nonlinear_shape_mode_error_decreases_on_refinement():
    ref=reference('broader')
    errs=[abs(rate('broader',bins=n)['Rphi']/ref['R_phi']-1)
          for n in (64,128,256,512)]
    assert all(a>b for a,b in zip(errs,errs[1:]))
    assert 3<errs[-2]/errs[-1]<5


def test_parent_momentum_tail_cutoff_is_independent_large_uncertainty():
    rows=cutoff_audit(cutoffs=(2.,4.,8.,12.,16.))
    vals=[abs(r['Rphi']) for r in rows]
    assert all(x>0 for x in vals)
    assert vals[1]/vals[0]>2
    assert abs(vals[-1]/vals[-2]-1)<.001
    assert all(r['entropy']>0 for r in rows)
