"""Conditional mass normalization across cap prescriptions in the source halo units."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import constants
from scipy.special import zeta
from halo_bose_primordial_cap_core_reference import integrate, NCRIT
from halo_bose_finite_lowered_mixing_reference import solve
from halo_collisionless_bose_covariance_cap_reference import minimum_cap


def mass_from_central(R,V,nu):
    units=np.sqrt(np.pi/2)*constants.hbar**3/(constants.G*nu*(R*1e3*constants.parsec)**2*(V*1e3))
    return units**.25*constants.c**2/constants.e


def run(path):
    d=json.loads(path.read_text())
    D=d['distance_Mpc'];err=d['distance_error_Mpc'];R0=.6;V0=44.1
    def mass(D,ratio):
        return float(mass_from_central(R0*D/49.,V0,zeta(1.5)*ratio))
    refs=[mass(v,1.) for v in (49.,D,D-err,D+err)]
    i0=d['inclination_deg']; ei=d['inclination_error_deg']
    nuisance=[]
    for dist in (D-err,D,D+err):
        for i in (i0-ei,i0,i0+ei):
            vel=V0*np.sin(np.deg2rad(i0))/np.sin(np.deg2rad(i))
            val=mass_from_central(R0*dist/49.,vel,zeta(1.5))
            nuisance.append(dict(D_Mpc=dist,inclination_deg=i,conditional_mass_eV=val))
    out=[]
    for F in (.5123174296,1.,100.,5000.,7000.,200000.):
        delta=np.log1p(1/F)
        capped=integrate(0.,delta)/NCRIT
        shifted=integrate(delta,delta)/NCRIT
        assert 0.<shifted<capped<1.
        entry=dict(cap=F,clipped_max_central_density_over_critical=capped,
                   smooth_shift_max_density_over_critical=shifted,
                   clipped_mass_normalization_at_catalog_D_eV=mass(D,capped),
                   smooth_thermal_mass_normalization_at_catalog_D_eV=mass(D,shifted),
                   clipped_mass_relative_to_critical=capped**(-.25),
                   shifted_mass_relative_to_critical=shifted**(-.25),
                   arbitrary_shape_fixed_critical_moments_relative_mass_floor=(2.387455588567837/F)**.25)
        if F==7000:
            arms=[]
            for W in (2.,5.,8.):
                _,_,n,_,_,_=solve(F,W,return_profile=True)
                r=n/NCRIT
                assert r<1
                arms.append(dict(W=W,lowered_n0_over_critical=r,
                                 corrected_mass_normalization_eV=mass(D,r)))
            entry['finite_lowered_controls']=arms
        out.append(entry)
    # The same normalization predicts the density at source scales.
    rho= (V0*1e3)**2/(4*np.pi*constants.G*(R0*1e3*constants.parsec)**2)
    vd=V0*np.sqrt(zeta(2.5)/zeta(1.5))
    rho_msun_pc3=rho*constants.parsec**3/1.98847e30
    mx=minimum_cap(rho_msun_pc3,[vd]*3,mass=refs[0]*constants.e/constants.c**2)
    assert abs(mx/2.387455588567837-1)<3e-13
    return dict(scope="Source angular scale and V0 held as a benchmark, over different explicitly selected "
                      "halo density prescriptions. No fitting, uncertainty propagation beyond direct "
                      "distance rescaling, particle measurement or self-consistent revised cosmology.",
                critical_mass_scale_eV=dict(D49=refs[0],catalog_D=refs[1],
                                           catalog_distance_minus_error=refs[2],
                                           catalog_distance_plus_error=refs[3]),
                inclination_and_distance_rescale_grid=nuisance,
                cap_comparisons=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('data',type=Path)
    print(json.dumps(run(p.parse_args().data),indent=2,allow_nan=False))
