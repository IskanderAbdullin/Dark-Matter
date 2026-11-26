"""Conditional *fixed-fermion* Yukawa scalar mode positivity and independent N/Q.

Does not test a self-consistent three-species evolution or cosmic relic.
"""
import math
import numpy as np
from scipy.optimize import brentq
from kinetics.phi_bose_scalar_fixed_bath_strong import (
    mode_coefficients,exact_mode_step,integrate_weak,
    firstcell_exact_delta_moments,threshold,
)
from kinetics.phi_bose_moment_realisability_boundary import (
    boundary_crossing,difference_moments,
)
from kinetics.phi_bose_physical_cutcell_weak import (
    analytic_weak,MPHI,shell,
)


def test_mode_equilibrium_exact_chemical_potential_and_positive_semigroup():
    for k in (.01,.05,.1,.4):
        S,D,S1,D1=mode_coefficients(k,2.,48,T1=1.,T2=1.,
                                    mu1=-.15,mu2=.05)
        be=1/math.expm1(k+.2)
        assert S>0 and D>0
        assert math.isclose(S/D,be,rel_tol=4e-14)
        assert math.isclose(S1/D1,be,rel_tol=4e-14)
        assert math.isclose(
            exact_mode_step(k,2.,1e5,be,48,T1=1.,T2=1.,
                            mu1=-.15,mu2=.05),be,rel_tol=3e-12)
    for k in (.005,.01,.05,.1,.4):
        b0=1/math.expm1(k)
        b1=exact_mode_step(k,2.,50.,b0,32)
        S,D,*_=mode_coefficients(k,2.,32)
        assert S>=0 and D>0 and b1>0 and b1<b0
        assert b1>=S/D
    # Fermi occupation inversion: D<0 allows positive Bose-stimulated growth.
    S,D,*_=mode_coefficients(.01,2.,32,T1=1.,T2=1.,mu1=.1,mu2=0.)
    assert S>0 and D<0
    assert exact_mode_step(.01,2.,1000.,1.,32,T1=1.,T2=1.,
                           mu1=.1,mu2=0.)>1.


def test_reversed_integration_independently_matches_all_N_and_Q():
    for P in (2.,8.,16.):
        reference=analytic_weak(P)
        errors=[]
        for n in (8,16,32):
            v=integrate_weak(P,16,n,n)
            error=max(abs(v[k]/reference[k]-1)
                      for k in ('Rphi','Qphi','Q1','Q2'))
            errors.append(error)
            assert abs(v['energy_defect'])<2e-24
            assert v['Rphi']<0 and v['Qphi']<0
        assert errors[2]<2e-9,(P,errors)
        assert errors[1]<errors[0] and errors[2]<errors[1]
    v=integrate_weak(2.,16,32,32)
    assert math.isclose(v['bycell'][0,0],
                        -2.762149360447387e-12,rel_tol=4e-11)
    assert math.isclose(v['bycell'][0,1],
                        -1.2277997619515968e-14,rel_tol=4e-11)


def test_kinematic_support_and_absence_of_artificial_zero_mode():
    low,hi=shell(2.)[1:3]
    assert threshold(low*.9,2.) is None
    assert threshold(hi*1.1,2.) is None
    assert threshold(low,2.)>=0
    assert threshold(.01,2.)>0
    for k in (low*.9,hi*1.1):
        assert mode_coefficients(k,2.)==(0.,0.,0.,0.)


def test_positive_actual_mode_evolution_is_outside_affine_eta_family():
    for n in (24,32,40):
        b=boundary_crossing(n)
        for fraction,sign in ((.8,-1),(1.1,+1),(2.,+1)):
            dt=b['boundary_time']*fraction
            moment_delta=firstcell_exact_delta_moments(n,2.,dt,48,32)
            def solve_right(er):
                return (difference_moments(b['upper'],0.,er)[1]
                        -moment_delta[1])
            er=brentq(solve_right,b['upper']*1e-5,
                      b['upper']+1e-5,xtol=1e-16)
            critical_left=difference_moments(b['upper'],0.,er)[0]
            gap=moment_delta[0]-critical_left
            assert sign*gap>1e-15,(n,fraction,gap)
            # Finite exact semigroup differs predictably from frozen-source
            # linearization; at 2t* the difference is below 0.1% here.
            assert abs(moment_delta[0]/(dt*b['C_left_right'][0])-1)<1e-3
            assert abs(moment_delta[1]/(dt*b['C_left_right'][1])-1)<1e-3
            for k in np.linspace(shell(2.)[1],b['upper'],9):
                if k<=MPHI:continue
                b0=1/math.expm1(k)
                assert exact_mode_step(float(k),2.,dt,b0,32)>0


def test_exact_mode_delta_derivative_matches_independent_initial_weak_source():
    for n in (24,32):
        b=boundary_crossing(n)
        for dt in (.01,1.,10.):
            delta=firstcell_exact_delta_moments(n,2.,dt,48,32)
            assert np.allclose(delta/dt,np.array(b['C_left_right']),
                               rtol=7e-5,atol=1e-25)
