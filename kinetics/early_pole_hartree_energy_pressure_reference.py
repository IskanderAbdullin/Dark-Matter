"""Conserving energy and pressure for the vacuum-subtracted Hartree comparison."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from early_pole_shell_hartree_feedback_reference import S,setup,pole


def moments(N,P,m0,lam,a):
    n=N/a**3
    p=P/a
    m=pole(n.sum(),lam,m0,p,n/n.sum())
    E=np.hypot(p,m)
    U=(m*m-m0*m0)**2/(2*lam) if lam else 0.
    K=n@E
    pr=n@(p*p/(3*E))+U
    rho=K-U
    assert rho>0 and pr>=0
    return rho,pr,m,U,K


def run(f):
    rows=json.loads(f.read_text())['examples']
    p,w,_=setup()
    out=[]
    for r in rows[1:5]:
        lam=r['lambda_self_quartic']
        x=r['completed_to_x']
        N=r['final_shell_yield']*S*w
        m0=1e-9*x
        base=moments(N,p,m0,lam,1.)
        assert abs(base[2]/r['final_scalar_pole_over_T']-1)<3e-10
        # Work-free diluted shell: verify the derivative of the conserving
        # rho against -3(rho+P), at every indicative stage.
        arc=[]
        for a in (1.,10.,1e3,1e6,1e7,1e8,1e9,1e10):
            old=moments(N,p,m0,lam,a*np.exp(-1e-4))
            z=moments(N,p,m0,lam,a)
            nex=moments(N,p,m0,lam,a*np.exp(1e-4))
            slope=(nex[0]-old[0])/2e-4
            dd=abs(slope/(-3*(z[0]+z[1]))-1)
            assert dd<1e-6
            v=(nex[1]-old[1])/(nex[0]-old[0])
            formula=(z[2]**2-m0*m0)/(3*z[2]**2-m0*m0)
            arc.append(dict(relative_scale=a,pressure_over_density=z[1]/z[0],
                            adiabatic_dP_dRho=v,zero_momentum_formula=formula,
                            collisionless_hold_particle_energy_fraction_removed_as_double_count=z[3]/z[4],
                            scalar_pole_over_vacuum_mass=z[2]/m0,
                            energy_conservation_derivative_defect=dd))
        val=moments(N,p,m0,lam,1.)[0]
        pn=N.copy()
        pn[3]*=1+1e-4
        mn=N.copy()
        mn[3]*=1-1e-4
        q=(moments(pn,p,m0,lam,1.)[0]-moments(mn,p,m0,lam,1.)[0])/(pn[3]-mn[3])
        assert abs(q/np.hypot(base[2],p[3])-1)<2e-7
        # At p=0 the gap is purely cubic; the condition m^2=2 m0^2
        # is mapped back to an illustrative expansion time while holding count.
        ae=(lam*N.sum()/(2*np.sqrt(2)*m0**3))**(1/3)
        assert abs(moments(N,p,m0,lam,ae)[2]/(np.sqrt(2)*m0)-1)<1e-4
        out.append(dict(lambda_self_quartic=lam,fermion_pole_over_T=r['fermion_pole_over_T'],
                        energy_source_functional_derivative_defect=abs(q/np.hypot(base[2],p[3])-1),
                        scale_factor_growth_to_self_mass_squared_equal_vacuum=ae,
                        frozen_examples=arc))
    return dict(scope="A variational Hartree consistency test for a held comoving shell and fixed vacuum "
                      "mass. Adiabatic background w and dP/dRho only; no perturbation closure, "
                      "no new decay history, other quantum corrections or observed constraint.",
                histories=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('shell_ledger',type=Path)
    print(json.dumps(run(p.parse_args().shell_ledger),indent=2,allow_nan=False))
