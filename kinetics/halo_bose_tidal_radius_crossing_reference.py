"""Initial energy-selected tidal loss corrected for angular-momentum crossing."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from scipy.stats import qmc
from halo_bose_finite_lowered_mixing_reference import solve
from halo_bose_finite_evaporation_kernel_reference import occupied
from halo_bose_orbital_elastic_clock_reference import R0,V0,SIGMA_OVER_M,constants,GYR


def local(y,d,W,xr,Z,angular_n):
    v=np.sqrt(2*(W-y))*np.r_[0.,np.geomspace(1e-7,1.,8000)]
    f=occupied(y,v*v/2,d,W)
    n=cumulative_simpson(v*v*f,x=v,initial=0.)
    good=np.r_[True,np.diff(n)>1e-16*n[-1]]
    q=n[good]/n[-1]
    finv=PchipInterpolator(q,v[good],extrapolate=False)
    v,u=finv(np.minimum(Z[:,:2],q[-1])).T
    mu=2*Z[:,2]-1
    S=v*v+u*u
    G=np.sqrt((u-v)**2+2*u*v*(1-mu))
    b=.25*np.sqrt(np.maximum(S*S-(2*u*v*mu)**2,0.))
    a=y+S/4
    e0=np.maximum(a+b-W,0.)
    cond=e0>0
    # Surplus threshold for *all* outgoing angular momenta to fit
    # through the outer radius, given the fixed scattering point.
    k=W-y
    xsq=xr*xr
    ecap=min(xsq*k/max(1-xsq,1e-14),2*W)
    M=np.minimum(ecap,e0[cond])
    ang,w=roots_legendre(angular_n)
    u=(1+ang)/2
    t=M[:,None]*(1-u*u)
    V=np.maximum(2*a[cond,None]-W-t-y,0.)
    # All directions of the high energy daughter are equally likely
    # relative to radial direction, after averaging the spherical bath.
    blocked=np.sqrt(np.maximum(0.,1.-t/(xsq*(k+t))))
    I=M*((1+occupied(y,V,d,W))*blocked@(u*w))/b[cond]
    return n[-1]/np.sqrt(2), G[cond]@I/len(G)


def control():
    p2=.4;k=2.;emax=1.1
    u,w=roots_legendre(56)
    M=min(p2*k/(1-p2),emax)
    ee=M*(1-((u+1)/2)**2)
    approx=M*(np.sqrt(1-ee/(p2*(k+ee)))@((u+1)*w/2))
    exact=quad(lambda e:np.sqrt(1-e/(p2*(k+e))),0.,M,epsabs=1e-13)[0]
    assert abs(approx/exact-1)<3e-12
    return abs(approx/exact-1)


def run(ledger):
    r0=V0*V0/(4*np.pi*constants.G*R0*R0)
    C=r0*SIGMA_OVER_M*V0*GYR
    checks=control()
    out=[]
    for row in json.loads(ledger.read_text())['cases']:
        F,W=row['peak_occupation'],row['escape_depth']
        d,h,n,N,sol,lr=solve(F,W,return_profile=True)
        extra=sol(lr)[1]
        arms=[]
        for A,P,S,L in ((66,14,641,16),(104,16,641,22),(104,16,734,22)):
            x,w=roots_legendre(A)
            z=-13+(lr+13)*(1+x)/2
            y=np.minimum(sol(z)[0],W-1e-13)
            Z=qmc.Sobol(3,scramble=True,seed=S).random_base2(P)
            T=np.array([local(a,d,W,np.exp(q-lr),Z,L) for a,q in zip(y,z,strict=True)])
            ds=w*np.exp(3*z)*(lr+13)/2
            assert abs((ds@T[:,0])/N-1)<8e-5
            arms.append(C*(ds*T[:,0]**2)@T[:,1]/(2*N*n))
        val=np.mean(arms[1:])
        dif=max(abs(arms[0]/val-1),abs(arms[1]/arms[2]-1))
        assert dif<.032,(W,arms,dif)
        J=row['tidal_flux_fraction_per_Gyr']
        I=row['infinity_flux_fraction_per_Gyr']
        assert I<J-val<J
        out.append(dict(peak_occupation=F,escape_depth=W,tidal_energy_selected_rate=J,
                        angular_momentum_blocked_from_crossing_radius_rate=val,
                        actual_outward_accessible_rate=J-val,
                        isolated_positive_energy_rate=I,
                        negative_energy_among_outward_accessible_fraction=(J-val-I)/(J-val),
                        grid_scramble_defect_of_blocked_rate=dif))
    return dict(scope="Initial elastic production in a fixed spherical lowered halo. "
                      "No rescattering en route, time-dependent tides or return dynamics. "
                      "Only reachability of the geometric outer sphere and the sign of isolated "
                      "orbital energy; not a time-dependent escape model.",
                analytic_selection_quadrature_defect=checks,outputs=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ledger',type=Path)
    print(json.dumps(run(p.parse_args().ledger),indent=2,allow_nan=False))
