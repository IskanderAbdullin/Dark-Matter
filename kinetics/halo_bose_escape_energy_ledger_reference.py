"""Initial mass and energy flux with tidal versus zero-at-infinity boundaries."""
import json
import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from scipy.stats import qmc
from halo_bose_finite_lowered_mixing_reference import solve
from halo_bose_finite_evaporation_kernel_reference import occupied
from halo_bose_orbital_elastic_clock_reference import R0,V0,SIGMA_OVER_M,constants,GYR


def polar(y,d,W,Vg,S,cut):
    b=Vg/2
    a=y+S/4
    C=np.clip((cut-a)/b,0.,1.)
    width=1-C
    base=np.maximum(a+b-cut,0.)
    # Energy above this particular cut: 0 at the lower polar endpoint.
    # Segmented quadrature isolates the residual narrow angular Bose peak.
    x,w=roots_legendre(8)
    z=np.r_[0.,np.geomspace(.000001,1.,10)]
    A,B=z[:-1,None],z[1:,None]
    u=((A+B)/2+(B-A)/2*x).ravel()
    w=((B-A)/2*w).ravel()
    Eout=a[:,None]+b[:,None]*(1.-width[:,None]*u)
    eL=np.maximum(2*a[:,None]-Eout-y,0.)
    fac=1.+occupied(y,eL,d,W)
    m=width*(fac@w)
    K=width*((fac*np.maximum(Eout-cut,0.))@w)
    # Independent identity for the vacuum angular moment.
    assert np.max(abs(width*base/2-width*(np.maximum(Eout-cut,0.)@w)))<2e-14
    return m,K


def local(y,d,W,Z,extra):
    g=np.r_[0.,np.geomspace(1e-7,1.,8000)]
    v=np.sqrt(2*(W-y))*g
    f=occupied(y,v*v/2,d,W)
    c=cumulative_simpson(v*v*f,x=v,initial=0.)
    keep=np.r_[True,np.diff(c)>1e-16*c[-1]]
    qq=c[keep]/c[-1]
    inv=PchipInterpolator(qq,v[keep],extrapolate=False)
    u,v=inv(np.minimum(Z[:,:2],qq[-1])).T
    mu=2*Z[:,2]-1
    S=u*u+v*v
    G=np.sqrt((u-v)**2+2*u*v*(1-mu))
    VG=.5*np.sqrt(np.maximum(S*S-(2*u*v*mu)**2,0.))
    vals=[]
    for cut in (W,W+extra):
        tot=np.zeros(2)
        for i in range(0,len(G),1536):
            k=slice(i,i+1536)
            use=(y+S[k]/4+VG[k]/2)>cut
            if not use.any():
                continue
            p,K=polar(y,d,W,VG[k][use],S[k][use],cut)
            gg=G[k][use]
            tot+=np.array([gg@p,gg@K])
        vals.extend(tot/len(G))
    assert vals[0]>=vals[2]-1e-12
    return np.r_[c[-1]/np.sqrt(2),vals]


def run():
    rho0=V0*V0/(4*np.pi*constants.G*R0**2)
    ans=[]
    for W in (2.,5.,8.):
        F=7000.
        d,h,n,N,sol,lr=solve(F,W,return_profile=True)
        extra=sol(lr)[1]
        arm=[]
        for num,p,seed in ((70,14,275),(104,16,275),(104,16,581)):
            x,w=roots_legendre(num)
            zz=-13+(1+x)*(lr+13)/2
            t=np.minimum(sol(zz)[0],W-1e-13)
            v=qmc.Sobol(3,scramble=True,seed=seed).random_base2(p)
            b=np.array([local(y,d,W,v,extra) for y in t])
            du=w*np.exp(3*zz)*(lr+13)/2
            nu=du@b[:,0]
            assert abs(nu/N-1)<8e-5
            coef=du*b[:,0]**2/(2*N*n)
            arm.append(coef@b[:,1:])
        lo,b,c=arm
        q=(b+c)/2
        defect=max(np.max(abs(lo/q-1)),np.max(abs(b/c-1)))
        fine=float(np.max(abs(b/c-1)))
        assert defect<.065 and fine<.012,(W,arm)
        J,K,J2,K2=q*rho0*SIGMA_OVER_M*V0*GYR
        E_t=K/J-extra
        E2=K2/J2
        ans.append(dict(peak_occupation=F,escape_depth=W,exterior_energy_offset=extra,
                        tidal_flux_fraction_per_Gyr=J,tidal_mean_escaping_energy_over_mV02_at_infinity_zero=E_t,
                        infinity_flux_fraction_per_Gyr=J2,infinity_mean_excess_energy_over_mV02=E2,
                        true_to_tidal_fraction_of_events=J2/J,
                        tidal_flux_of_particle_energy_per_NmV02_per_Gyr=J*E_t,
                        isolated_flux_of_particle_energy_per_NmV02_per_Gyr=K2,
                        maximum_grid_and_scramble_defect=defect,
                        fine_scramble_defect=fine))
    return dict(scope="Fixed isolated spherical vacuum exterior for the infinity energy offset; "
                      "tidal removal is a distinct imposed sink. The exported particle energy is "
                      "an instantaneous flux diagnostic in the held potential, not a recalculated "
                      "self-gravitational energy, orbit return, or evolution law.",cases=ans)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
