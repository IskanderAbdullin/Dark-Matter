"""Sharp local density bound at finite occupation and fixed velocity covariance."""
import json
import numpy as np
from scipy import constants
from scipy.integrate import quad
from scipy.special import zeta

MSUN=1.98847e30
M=86*constants.e/constants.c**2

def minimum_cap(rho_msun_pc3,sigmas_km_s,mass=M):
    rho=rho_msun_pc3*MSUN/constants.parsec**3
    prod=np.prod(np.asarray(sigmas_km_s)*1e3)
    return rho*6*np.pi**2*constants.hbar**3/(5**1.5*mass**4*prod)


def run():
    # In transformed coordinates with unit covariance, a filled
    # ball of radius sqrt(5) exactly saturates the bound.
    z=np.sqrt(5.)
    vol=quad(lambda x:4*np.pi*x*x,0.,z)[0]
    norm=quad(lambda x:4*np.pi*x**4,0.,z)[0]/vol
    assert abs(norm-3)<1e-13
    assert abs(vol/(4*np.pi*5**1.5/3)-1)<1e-13
    references=[]
    for rho in (.01,.1,1.):
        for vel in (10.,20.,50.):
            references.append(dict(scalar_density_Msun_per_pc3=rho,sigmas_km_s=[vel]*3,
                                   minimum_required_occupation_cap=minimum_cap(rho,[vel]*3)))
    aniso=[]
    for ss in ((20.,20.,20.),(20.,20.,10.),(30.,20.,10.),(30.,30.,10.)):
        aniso.append(dict(scalar_density_Msun_per_pc3=.1,sigmas_km_s=ss,
                          minimum_required_occupation_cap=minimum_cap(.1,ss)))
    const=(6*np.pi**2/(5**1.5*(2*np.pi)**1.5))*(zeta(1.5)**2.5/zeta(2.5)**1.5)
    assert 1.<const<10.
    return dict(scope="Local collisionless necessary phase-space moment bound, one scalar degree of "
                      "freedom, Newtonian velocities, no observed halo density or velocity inference. "
                      "Not sufficient for a halo phase-space profile or dynamical history.",
                mass_eV_over_c2=86.,
                isotropic_fiducial_examples=references,anisotropic_fiducial_examples=aniso,
                critical_ideal_bose_moments_only_minimum_cap=const)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
