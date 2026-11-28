"""Local one-mode linearized Bose elastic shear dissipation, not an eigenvalue."""
import json
import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma, ndtri, roots_legendre
from scipy.stats import qmc
from halo_bose_quantum_orbital_elastic_reference import average
from halo_bose_orbital_elastic_clock_reference import (
    background, R0, V0, SIGMA_OVER_M, constants, GYR,
)


def radial(a0):
    v=np.r_[0.,np.geomspace(1e-7,14.,12500)]
    z=a0+v*v/2
    f=np.exp(-z)/(-np.expm1(-z))
    cum=cumulative_simpson(v*v*f,x=v,initial=0.)
    keep=np.r_[True,np.diff(cum)>1e-16*cum[-1]]
    x=cum[keep]/cum[-1]
    inv=PchipInterpolator(x,v[keep],extrapolate=False)
    J=lambda fun:quad(lambda u:fun(u,np.exp(-a0-u*u/2)/(-np.expm1(-a0-u*u/2))),
                     0.,14.,epsabs=1e-11,epsrel=1e-9)[0]
    q=J(lambda u,f:u*u*f)
    assert abs(q/cum[-1]-1)<2e-8
    chi=J(lambda u,f:u**6*f*(1+f))/(15*q)
    var=J(lambda u,f:u**4*f)/(3*q)
    return lambda x:inv(np.minimum(x,max(cum[keep])/cum[-1])),q,chi,var


def angles(n):
    x,w=roots_legendre(n)
    # Resolve endpoint moments where the outgoing energies are far apart.
    edge=1-np.r_[np.geomspace(1.,1e-5,10),0.]
    lo,hi=edge[:-1,None],edge[1:,None]
    return (((lo+hi)/2+(hi-lo)/2*x).ravel(),((hi-lo)/2*w).ravel())


def qmc_bracket(a0,rad,seed,power,order):
    r=qmc.Sobol(3,scramble=True,seed=seed).random_base2(power)
    uv=rad(r[:,:2]);v,u=uv[:,0],uv[:,1]
    x=2*r[:,2]-1
    S=v*v+u*u
    g=np.sqrt((u-v)**2+2*u*v*(1-x))
    V=.5*np.sqrt((u-v)**2+2*u*v*(1+x))
    alpha2=np.minimum(((v*v-u*u)/(2*V*g))**2,1.)
    b=V*g/2
    a=a0+S/4
    xmin=a0+(v*u*x)**2/(S+4*b)
    c,w=angles(order)
    J=[]
    for i in range(0,len(a),768):
        sl=slice(i,i+768)
        z=xmin[sl,None]+b[sl,None]*(1-c)
        W=1/((-np.expm1(-z))*(-np.expm1(-2*a[sl,None]+z)))
        L=W@(w*(1-c*c))
        M=average(a[sl],b[sl],xmin[sl])
        res=(1.-alpha2[sl])*M+.5*(3*alpha2[sl]-1)*L
        # This is an integral of 1-cos(angle_initial_to_final)^2.
        assert np.min(res)>0
        J.extend(g[sl]**5*res)
    return np.mean(J)/80.


def geometric_test():
    # Alternative orientations: condition on random cosines and independently
    # reconstruct full Cartesian outgoing and incoming in the center-of-mass.
    rng=np.random.default_rng(17512)
    n=45000
    z=rng.normal(size=(n,3));z/=np.linalg.norm(z,axis=1)[:,None]
    q=rng.normal(size=(n,3));q/=np.linalg.norm(q,axis=1)[:,None]
    a=np.asarray([.1,.22,-.2])
    g1=z*1.223;g2=q*1.223
    v=a+g1/2;u=a-g1/2;v2=a+g2/2;u2=a-g2/2
    assert np.max(abs(np.sum(v*v+u*u-v2*v2-u2*u2,axis=1)))<1e-14
    psi=lambda v:v[:,0]*v[:,1]
    actual=np.mean((psi(v)+psi(u)-psi(v2)-psi(u2))**2)
    expect=1.223**4/30.
    assert abs(actual/expect-1)<.018
    return actual/expect


def run():
    rel=geometric_test()
    rho0=V0**2/(4*np.pi*constants.G*R0**2)
    clock=rho0*SIGMA_OVER_M*V0*GYR
    rows=[]
    for F in (100.,7000.):
        d=np.log1p(1/F)
        pot,n0=background(d)
        for r in (.1,.3,1.,3.,10.,20.):
            y=pot(r)[0]
            fn,n,chi,var=radial(y+d)
            lo=np.array([qmc_bracket(y+d,fn,k,15,8) for k in (17,39,51)])
            hi=np.array([qmc_bracket(y+d,fn,k,17,12) for k in (17,39,51)])
            err=max(np.ptp(hi),abs(lo.mean()-hi.mean()))/hi.mean()
            assert err<.018, (F,r,lo,hi)
            lam=hi.mean()/chi
            nc=n/(np.sqrt(2)*n0)
            mv=8/(5*gamma(1.5))*np.sqrt(var)
            rows.append(dict(central_occupation=F,radius_R=r,local_density_over_central=nc,
                             shear_Rayleigh_per_Gyr=lam*clock*nc,
                             matched_Maxwell_vacuum_Rayleigh_per_Gyr=mv*clock*nc,
                             shear_over_matched_vacuum_Maxwell=lam/mv,
                             combined_scramble_and_resolution_defect=err))
    return dict(scope="Linearized local ideal Uehling-Uhlenbeck shear test function only. "
                      "Matched Maxwell comparison uses the same local number and velocity variance. "
                      "No spectral gap, actual mode, nonlinear transport, or global equilibration time.",
                independent_cartesian_shear_tensor_test_ratio=rel,examples=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
