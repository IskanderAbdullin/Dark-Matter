"""Local elastic dissipation of an isotropic nonthermal fourth-moment distortion."""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
from scipy.stats import qmc
from halo_bose_linearized_shear_reference import angles,radial
from halo_bose_orbital_elastic_clock_reference import (
    background, R0, V0, SIGMA_OVER_M, constants, GYR,
)


def susceptibility(a):
    def integrate(j):
        def g(u):
            exp=np.exp(-a-u*u/2)
            f=exp/(-np.expm1(-a-u*u/2))
            return u**(2+2*j)*f*(1+f)
        return quad(g,0.,14.,epsabs=1e-10,epsrel=3e-10)[0]
    H=[integrate(j) for j in range(5)]
    coeff=np.linalg.solve([[H[0],H[1]],[H[1],H[2]]],[H[2],H[3]])
    C0,C2=coeff
    assert abs(H[2]-C0*H[0]-C2*H[1])<2e-11*H[2]
    assert abs(H[3]-C0*H[1]-C2*H[2])<2e-11*H[3]
    val=H[4]-C0*H[2]-C2*H[3]
    assert val>0.
    return val,C0,C2


def bracket(a0,fn,seed,L,M):
    x=qmc.Sobol(3,scramble=True,seed=seed).random_base2(L)
    v,u=fn(x[:,:2]).T
    mu=2*x[:,2]-1
    S=v*v+u*u
    g=np.sqrt((u-v)**2+2*u*v*(1-mu))
    V=.5*np.sqrt((u-v)**2+2*u*v*(1+mu))
    b=V*g/2
    a=a0+S/4
    xmin=a0+(u*v*mu)**2/(S+4*b)
    al=np.minimum(((v*v-u*u)/(2*V*g))**2,1.)
    c,w=angles(M)
    J=0.
    for i in range(0,len(a),1024):
        s=slice(i,i+1024)
        z=xmin[s,None]+b[s,None]*(1-c)
        W=1/((-np.expm1(-z))*(-np.expm1(-2*a[s,None]+z)))
        test=(((c*c-al[s,None])**2)*W)@w
        J+=np.dot(g[s]*(V[s]*g[s])**4,test)
    return J/len(a)


def run():
    rho0=V0**2/(4*np.pi*constants.G*R0**2)
    clock=rho0*SIGMA_OVER_M*V0*GYR
    out=[]
    for F in (100.,7000.):
        d=np.log1p(1/F)
        pot,n0=background(d)
        for r in (.1,.3,1.,3.,10.,20.):
            z=pot(r)[0]+d
            fn,n,_,var=radial(z)
            chi,c0,c2=susceptibility(z)
            lo=np.array([bracket(z,fn,s,15,8) for s in (17,39,51)])
            hi=np.array([bracket(z,fn,s,17,12) for s in (17,39,51)])
            err=max(np.ptp(hi),abs(lo.mean()-hi.mean()))/hi.mean()
            assert err<.025,(F,r,lo,hi)
            lam=hi.mean()/(chi/n)
            vac=16/(15*gamma(1.5))*np.sqrt(var)
            out.append(dict(central_occupation=F,radius_R=r,coeff_0=c0,coeff_v2=c2,
                            shape_Rayleigh_per_Gyr=lam*clock*n/(np.sqrt(2)*n0),
                            matched_Maxwell_shape_Rayleigh_per_Gyr=vac*clock*n/(np.sqrt(2)*n0),
                            shape_over_matched_Maxwell=lam/vac,
                            combined_scramble_and_resolution_defect=err))
    # Exact Maxwell invariant projection in unit one-dimensional variance.
    e=[1.,3.,15.,105.,945.]
    c0,c2=np.linalg.solve([[e[0],e[1]],[e[1],e[2]]],[e[2],e[3]])
    assert abs(e[4]-c0*e[2]-c2*e[3]-120)<1e-12
    assert abs(c0+15)<1e-12 and abs(c2-10)<1e-12
    return dict(scope="Local scalar-energy-conserving fourth-moment response of the ideal quantum "
                      "elastic operator; no global mode, equilibration proof or cosmological evolution.",
                maxwell_projected_denominator_unit_variance=120.,rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
