"""Regression tests for independent momentum grids and parent-angle kernel."""
import numpy as np
from scipy.special import expit
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'kinetics'))
from resolved_parent_bath import ResolvedParent
from parent_angle_kernel import ParentAngle,positive_hat_integrals

def test_hat_basis_number_and_first_moment():
    nodes=np.array([.1,.2,.5,1.,3.])
    a=np.array([.11,.2,.4,1.1]);b=np.array([.4,.5,2.,2.8])
    w=positive_hat_integrals(nodes,a,b)
    assert np.min(w)>=0
    assert np.allclose(w.sum(axis=0),b-a,rtol=1e-14)
    assert np.allclose(nodes@w,(b*b-a*a)/2,rtol=1e-14)

def test_parent_grid_aliasing_counterexample():
    kw=dict(n=120,nparent=120,kappa=.006,scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8)
    values=[]
    for grid in ['legacy','uniform']:
        m=ResolvedParent(parent_grid=grid,**kw);m.geometry(.2)
        load=np.bincount(np.r_[m.lo,m.hi],weights=np.r_[m.wl*m.event_rate0,m.wh*m.event_rate0],minlength=m.n1)
        values.append(np.max(load/(m.v1*m.kappa*.2*m.t1)))
    assert values[0]>50 and values[1]<5,values

def test_angle_detailed_balance_positivity_and_normalization():
    m=ParentAngle(n=80,nparent=80,kappa=.006,scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8,
                  parent_grid='balanced',fermion_grid='balanced',parent_floor=0.)
    r=m.validate()
    assert r['balance']<1e-10 and r['energy']<1e-12
    assert r['max_shell_excess']<1e-10
    assert r['max_projection_energy_error']<1e-12
    m.geometry(.2)
    y=np.r_[np.full(m.n,.2),expit(-m.energy_chi),expit(-m.e1)]
    c=m.dark(y);assert np.isfinite(c).all()

def test_angle_regularization_small_enough_for_late_detailed_balance():
    kw=dict(n=80,nparent=80,kappa=.006,scalar_ratio=1e-12,fermion_ratio=.9,qmin=1e-8,
            parent_grid='balanced',fermion_grid='balanced')
    vals=[]
    for floor in [1e-12,1e-18]:
        m=ParentAngle(parent_floor=floor,**kw);m.geometry(20.)
        eq=expit(-m.energies);eq[:m.n]=1/np.expm1(m.energy_phi)
        vals.append(float(abs(m.dark(eq))@m.vols/(m.event_rate0@expit(-m.e1[m.jp]))))
    assert vals[0]>1e-5 and vals[1]<1e-8,vals


def test_separate_rectangle_same_grid_closure():
    from parent_angle_kernel import SeparatedRectangle
    from separated_grid_handoff import GridClosure
    m=SeparatedRectangle(n=48,nparent=48,kappa=.00685,scalar_ratio=1e-12,
                         fermion_ratio=.9,parent_grid='uniform',fermion_grid='balanced',
                         parent_floor=0.,qmin=1e-8)
    m.geometry(.4)
    eq=expit(-m.energies);eq[:m.n]=1/np.expm1(m.energy_phi)
    err=np.linalg.norm(m.dark(eq)*m.vols,ord=1)
    assert err<1e-10
    assert GridClosure(m).validate()['jacobian_relative_error']<1e-5
