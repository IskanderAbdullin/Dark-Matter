"""Physical Yukawa continuum mode bounds and independent nonthermal weak ledgers."""
import math
import numpy as np
from scipy.special import expit
from kinetics.phi_bose_coupled_continuum_modes import (
    fd1,fd2,bephi,parent_mode,daughter_mode,scalar_mode,
    integrated_parent_weak,integrated_scalar_weak,
)
from kinetics.phi_bose_physical_cutcell_weak import M1


def constant(n):
    return lambda E: np.zeros_like(np.asarray(E),dtype=float)+n


def test_inward_fermion_and_bose_boundaries_nonthermal():
    for p in (.0,.1,1.,4.):
        gain,loss=parent_mode(p,fd2,bephi,36)
        assert gain>=0 and loss>=0
        assert gain*(1-0)-loss*0>=0
        assert gain*(1-1)-loss*1<=0
    for P in (2.,8.,16.):
        for e2 in (.03,.052,.1,.4):
            gain,loss=daughter_mode(e2,P,fd1,bephi,36)
            assert gain>=0 and loss>=0
            assert gain*(1-0)-loss*0>=0
            assert gain*(1-1)-loss*1<=0
        for k in (.002,.01,.05,.2):
            S,D=scalar_mode(k,P,fd1,fd2,36)
            assert S>=0
            assert S-D*0>=0
    assert parent_mode(.5,constant(0.),constant(0.))[0]==0
    S,D=scalar_mode(.05,2.,constant(0.),constant(.6))
    assert S==0 and D>0
    S,D=scalar_mode(.05,2.,constant(.9),constant(.1))
    assert S>0 and D<0


def test_shared_nonthermal_parent_N_three_Q_H():
    for P in (2.,8.,16.):
        result=integrated_parent_weak(P,fd1,fd2,bephi,32,32)
        assert result['minW']>0
        assert result['Rphi']<0
        assert result['Qphi']<0
        assert result['Q1']>0
        assert result['Q2']<0
        assert abs(result['energy_defect'])<1e-24
        assert result['entropy']>0


def test_independent_scalar_first_nonthermal_convergence():
    for P in (2.,8.,16.):
        parent=integrated_parent_weak(P,fd1,fd2,bephi,96,96)
        errors=[]
        for n in (16,32,48,64,96):
            other=integrated_scalar_weak(P,fd1,fd2,bephi,n,n)
            errors.append(max(abs(other[0]/parent['Rphi']-1),
                              abs(other[1]/parent['Qphi']-1)))
        assert errors[-1]<1e-8,(P,errors)
        assert errors[-1]<errors[0],(P,errors)
        if P==16:
            assert errors[0]>1e-3


def test_common_T_chemical_equilibrium_mode_by_mode():
    T=1.;mu1=-.15;mu2=.05;muphi=-.2
    f1=lambda E: expit(-(np.asarray(E)-mu1)/T)
    f2=lambda E: expit(-(np.asarray(E)-mu2)/T)
    b=lambda E: 1/np.expm1((np.asarray(E)-muphi)/T)
    for P in (2.,8.,16.):
        for p in (.01,.1,1.,2.):
            G,L=parent_mode(p,f2,b,32)
            target=f1(math.hypot(p,M1))
            assert math.isclose(G/(G+L),float(target),rel_tol=3e-12)
        for e2 in (.025,.052,.1,.4):
            G,L=daughter_mode(e2,P,f1,b,40)
            assert G>0 and L>0
            assert math.isclose(G/(G+L),float(f2(e2)),rel_tol=3e-12)
        for k in (.002,.01,.05,.2):
            S,D=scalar_mode(k,P,f1,f2,40)
            assert S>=0 and D>0
            assert math.isclose(S/D,float(b(k)),rel_tol=3e-12)
        totals=integrated_parent_weak(P,f1,f2,b,32,32)
        assert abs(totals['Rphi'])<2e-25
        assert abs(totals['Qphi'])<2e-25
        assert abs(totals['Q1'])<2e-25
        assert abs(totals['Q2'])<2e-25
        assert abs(totals['entropy'])<1e-30
