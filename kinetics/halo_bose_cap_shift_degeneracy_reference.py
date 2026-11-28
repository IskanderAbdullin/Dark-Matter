"""Matched central-density cap-versus-shift controls in spherical Bose cores."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from halo_bose_primordial_cap_core_reference import integrate


def run():
    out=[]
    radius=np.array([.01,.03,.1,.3,1.,2.])
    for F in (100.,5000.,7000.,200000.):
        d=np.log1p(1/F)
        den0=integrate(0.,d)
        z=brentq(lambda x:integrate(x,x)-den0,d/100,d,xtol=1e-15)
        # Correct numerical cutoff in shifted evaluation: y+z always >= z.
        assert abs(integrate(z,z)/den0-1)<8e-12
        m0=-.5*integrate(0.,d,-.5)
        m1=-.5*integrate(z,z,-.5)
        dgrid=np.r_[0.,np.geomspace(d*1e-6,10.,900),d]
        dgrid=np.unique(dgrid)
        fa=PchipInterpolator(dgrid,[integrate(y,d)/den0 for y in dgrid])
        fb=PchipInterpolator(dgrid,[integrate(y+z,z)/den0 for y in dgrid])
        arms=[]
        for b0 in (0.,1.,5.):
            bs=.3
            bary=lambda x:b0/(1+(x/bs)**2)**2.5
            if b0:
                mass_b=lambda x:b0*x**3/(3*(1+(x/bs)**2)**1.5)
            else:
                mass_b=lambda x:np.zeros_like(x)
            st=1e-6
            def solve(f):
                zz=solve_ivp(lambda x,u:[u[1]/x**2,x*x*(f(u[0])+bary(x))],
                             (st,2.1),[(1+b0)*st*st/6,(1+b0)*st**3/3],
                             rtol=4e-10,atol=1e-13,dense_output=True)
                assert zz.success
                U,M=zz.sol(radius)
                return U,M
            ya,ma=solve(fa);yb,mb=solve(fb)
            dens=fa(ya)/fb(yb)-1
            total=np.sqrt(ma/mb)-1
            dark=np.sqrt((ma-mass_b(radius))/(mb-mass_b(radius)))-1
            check=None
            if F==7000.:
                _,ma0=solve(lambda y:integrate(y,d)/den0)
                _,mb0=solve(lambda y:integrate(y+z,z)/den0)
                check=float(np.max(abs(total-(np.sqrt(ma0/mb0)-1))))
                assert check<3e-8
            arms.append(dict(central_baryon_density_over_scalar=b0,baryon_scale_R=bs,
                             maximum_independent_direct_quadrature_speed_ratio_defect=check,
                             radius_R=radius.tolist(),clipped_over_shifted_scalar_density_minus_one=dens.tolist(),
                             clipped_over_shifted_total_speed_minus_one=total.tolist(),
                             clipped_over_shifted_scalar_speed_minus_one=dark.tolist()))
        out.append(dict(cap=F,cap_exponent=d,shifted_exponent_matching_central_density=z,
                        shifted_to_cap_exponent=z/d,shifted_peak_occupation_over_original_cap=1/(F*np.expm1(z)),
                        shifted_to_clipped_central_density_gradient_in_potential=m1/m0,
                        arms=arms))
    assert abs(out[-1]['shifted_to_cap_exponent']/((4/(3*np.pi))**2)-1)<2e-3
    assert abs(out[-1]['shifted_to_clipped_central_density_gradient_in_potential']/(3*np.pi*np.pi/16)-1)<2e-3
    return dict(scope="Matched spherical local toy halos at fixed actual scalar n0 and velocity scale; "
                      "some include a static Plummer density. No observation or quantitative "
                      "statement about line-of-sight gas, parameters, or actual profile inference.",
                high_cap_asymptote_shifted_to_cap_exponent=(4/(3*np.pi))**2,
                high_cap_asymptote_shifted_to_clipped_central_gradient=3*np.pi*np.pi/16,
                cases=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
