"""Entropy production from internal bound-bound elastic collisions of a lowered halo."""
import json
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
from scipy.stats import qmc
from halo_bose_finite_lowered_mixing_reference import solve
from halo_bose_finite_evaporation_kernel_reference import occupied
from halo_bose_orbital_elastic_clock_reference import R0,V0,SIGMA_OVER_M,constants,GYR


def local(y,d,W,r):
    cap=np.sqrt(2*(W-y))
    v=cap*np.r_[0.,np.geomspace(1e-7,1.,8500)]
    f=occupied(y,v*v/2,d,W)
    c=cumulative_simpson(v*v*f,x=v,initial=0.)
    keep=np.r_[True,np.diff(c)>1e-16*c[-1]]
    k=c[keep]/c[-1]
    inv=PchipInterpolator(k,v[keep],extrapolate=False)
    u,vv=inv(np.minimum(r[:,:2],k[-1])).T
    mu=2*r[:,2]-1
    S=u*u+vv*vv
    G=np.sqrt((u-vv)**2+2*u*vv*(1-mu))
    b=.25*np.sqrt(np.maximum(S*S-(2*u*vv*mu)**2,0.))
    k0=S/4
    capc=np.clip((W-y-k0)/b,0.,1.)
    cp=capc*r[:,3]
    f1=occupied(y,u*u/2,d,W)
    f2=occupied(y,vv*vv/2,d,W)
    f3=occupied(y,k0+b*cp,d,W)
    f4=occupied(y,np.maximum(k0-b*cp,0.),d,W)
    use=capc>0
    factors=lambda f:np.log(f[use])-np.log1p(f[use])
    A=factors(f1)+factors(f2)-factors(f3)-factors(f4)
    Wt=G[use]*capc[use]*(1+f3[use])*(1+f4[use])
    # Multiplicative bounded symmetrization: (L-G)/(L+G) = tanh(A/2).
    # Since A*tanh(A/2) is unchanged under reversing an event,
    # integrating L+G times that even factor equals twice integrating L.
    rate1=Wt@(A*np.tanh(A/2))/(2*len(r))
    rate2=Wt@A/(2*len(r))
    n=c[-1]/np.sqrt(2)
    return n,rate1,rate2


def run():
    rho0=V0*V0/(4*np.pi*constants.G*R0**2)
    clock=rho0*SIGMA_OVER_M*V0*GYR
    results=[]
    for W in (2.,5.,8.):
        F=7000.
        d,h,n,N,sol,lr=solve(F,W,return_profile=True)
        vals=[]
        for Nrad,p,seed in ((64,15,819),(104,17,819),(104,17,844)):
            x,w=roots_legendre(Nrad)
            l=-13+(lr+13)*(1+x)/2
            y=np.minimum(sol(l)[0],W-1e-12)
            r=qmc.Sobol(4,scramble=True,seed=seed).random_base2(p)
            arr=np.array([local(t,d,W,r) for t in y])
            rw=(lr+13)/2*w*np.exp(3*l)
            assert abs((rw@arr[:,0])/N-1)<8e-5
            vals.append((rw*arr[:,0]**2/(N*n))@arr[:,1:]*clock)
        g=np.mean(vals[1:],axis=0)
        assert g[0]>0 and abs(g[1]/g[0]-1)<.065,(W,vals)
        err=max(abs(vals[0][0]/g[0]-1),abs(vals[1][0]/vals[2][0]-1))
        assert err<.05,(W,vals,err)
        results.append(dict(peak_occupation=F,escape_depth=W,
                            initial_internal_entropy_production_per_particle_per_Gyr=g[0],
                            symmetrized_independent_form_check=g[1],
                            grid_and_scramble_defect=err))
    return dict(scope="Fixed lowered Bose halo. Local quantum elastic collisions retained only if "
                      "both incoming and both outgoing particles stay below the imposed tidal energy. "
                      "This sub-kernel conserves number and energy; includes no loss entropy, changing "
                      "phase distribution, spatial streaming in the operator or condensate dynamics.",
                outputs=results)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
