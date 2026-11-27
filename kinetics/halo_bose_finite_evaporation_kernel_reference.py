"""Direct two-body evaporation clock for fixed finite lowered Bose phase space."""
import json
import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from scipy.stats import qmc
from halo_bose_finite_lowered_mixing_reference import solve
from halo_bose_orbital_elastic_clock_reference import R0,V0,SIGMA_OVER_M,constants,GYR


def occupied(y,t,d,W):
    e=y+t
    return (np.exp(d+e)*np.expm1(np.maximum(W-e,0.))/
            (np.expm1(d+W)*np.expm1(d+e)))


def local(y,d,W,r):
    # Local independent pairs sampled from the exact held lowered initial DF.
    v=np.sqrt(2*(W-y))*np.r_[0.,np.geomspace(1e-7,1.,8000)]
    f=occupied(y,v*v/2,d,W)
    q=cumulative_simpson(v*v*f,x=v,initial=0.)
    nn=q[-1]
    good=np.r_[True,np.diff(q)>1e-16*nn]
    z=q[good]/nn
    fn=PchipInterpolator(z,v[good],extrapolate=False)
    v,u=fn(np.minimum(r[:,:2],z[-1])).T
    mu=2*r[:,2]-1
    S=v*v+u*u
    g=np.sqrt((u-v)**2+2*u*v*(1-mu))
    V=.5*np.sqrt((u-v)**2+2*u*v*(1+mu))
    b=V*g/2
    a=d+y+S/4
    C=np.minimum(np.maximum((W-y-S/4)/b,0.),1.)
    width=1-C
    # Kinetic energy of the remaining bound partner, after arbitrary sign
    # of outgoing label. Compute exponent's minimum without cancellation.
    low=d+y+(u*v*mu)**2/(S+4*b)
    inc=b*width
    integral=np.log1p(np.exp(-low)*(-np.expm1(-inc))/(-np.expm1(-low)))/b
    bose=width*(1-1/np.expm1(d+W))+integral
    assert np.min(bose-width)>-1e-14
    return nn/np.sqrt(2),float(g@width/len(r)),float(g@bose/len(r))


def polar_control():
    e=.8;d=.015;W=2.;y=.2;b=.55
    c=.43
    lo=d+y+e-b
    h=b*(1-c)
    exact=(1-c)*(1-1/np.expm1(d+W))+np.log1p(
        np.exp(-lo)*(-np.expm1(-h))/(-np.expm1(-lo)))/b
    control=quad(lambda t:1+occupied(y,e-b*t,d,W),c,1.)[0]
    assert abs(exact/control-1)<1e-14
    return exact/control-1


def run():
    test=polar_control()
    rho0=V0**2/(4*np.pi*constants.G*R0**2)
    clock=rho0*SIGMA_OVER_M*V0*GYR
    universe=[]
    for F,W in ((100.,2.),(7000.,2.),(7000.,5.),(7000.,8.)):
        d,h,n,N,solution,lrt=solve(F,W,return_profile=True)
        tests=[]
        for nr,bits,seed in ((82,14,194),(124,16,194),(124,16,177)):
            # Jacobi-free global radial integral in ln(r/R); raw mass check is independent.
            x,w=roots_legendre(nr)
            ll=-13+(x+1)*(lrt+13)/2
            y=np.minimum(solution(ll)[0],W-1e-13)
            S=qmc.Sobol(3,scramble=True,seed=seed).random_base2(bits)
            out=np.array([local(t,d,W,S) for t in y])
            xi2dx=w*np.exp(3*ll)*(lrt+13)/2
            mass=xi2dx@out[:,0]
            assert abs(mass/N-1)<8e-5
            measure=xi2dx*out[:,0]**2/(2*N*n)
            tests.append(np.r_[measure@out[:,1:],mass/N])
        a,b,c=tests
        dif=max(np.max(abs(a[:2]/b[:2]-1)),np.max(abs(c[:2]/b[:2]-1)))
        assert dif<.03,(F,W,tests,dif)
        rate=(b[:2]+c[:2])/2*clock
        universe.append(dict(central_occupation=F,escape_depth=W,tidal_radius_R=np.exp(lrt),
                             immediate_vacuum_fractional_escape_per_Gyr=rate[0],
                             with_bound_daughter_Bose_factor_per_Gyr=rate[1],
                             stimulated_to_vacuum_escape_ratio=rate[1]/rate[0],
                             joint_radial_and_qmc_relative_difference=dif,
                             mass_closure_error=max(abs(t[2]-1) for t in tests)))
    return dict(scope="Exact local elastic ejection of one particle per collision from a held stationary "
                      "lowered finite halo, with no returning debris. Fast initial sink diagnostics only; "
                      "no self-consistent distribution evolution, evaporation lifetime, or measured cross section.",
                polar_control_relative_defect=test,outputs=universe)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
