"""Ballistic duration outside a finite halo threshold in an isolated Kepler exterior."""
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

TARGETS=(1.,10.,100.)


def duration(r_t,D,eps,L2):
    mu=r_t*D
    A=mu/(2*eps)
    e=np.sqrt(np.maximum(1.-2*eps*L2/mu**2,1.e-30))
    x=np.clip((1-r_t/A)/e,-1.,1.)
    # Evaluate pi-anomaly without subtracting nearly equal angles.
    eta_c=np.arccos(-x)
    return 2*np.sqrt(A**3/mu)*(eta_c+e*np.sqrt(np.maximum(1-x*x,0.)))*R0/V0/GYR


def control():
    rt=5.;extra=.53;ep=.24;L2=.6
    t=duration(rt,extra,ep,L2)
    mu=rt*extra
    A=mu/(2*ep)
    rmax=A*(1+np.sqrt(1-2*ep*L2/mu**2))
    assert rmax>rt
    # Direct integration r=rmax-(rmax-rt)x^2 removes outer turning singularity.
    vr=lambda r:np.sqrt(2*(mu/r-ep)-L2/r**2)
    q=quad(lambda x:4*(rmax-rt)*x/vr(rmax-(rmax-rt)*x*x),0.,1.,epsabs=1e-10)[0]*R0/V0/GYR
    assert abs(t/q-1)<2e-10
    return abs(t/q-1)


def local(y,d,W,xr,rt,extra,Z):
    v=np.sqrt(2*(W-y))*np.r_[0.,np.geomspace(1e-7,1.,8000)]
    f=occupied(y,v*v/2,d,W)
    n=cumulative_simpson(v*v*f,x=v,initial=0.)
    good=np.r_[True,np.diff(n)>1e-16*n[-1]]
    q=n[good]/n[-1]
    fn=PchipInterpolator(q,v[good],extrapolate=False)
    v,u=fn(np.minimum(Z[:,:2],q[-1])).T
    cos=2*Z[:,2]-1
    S=v*v+u*u
    G=np.sqrt((u-v)**2+2*u*v*(1-cos))
    b=.25*np.sqrt(np.maximum(S*S-(2*u*v*cos)**2,0.))
    a=y+S/4
    C=np.clip((W-a)/b,0.,1.)
    hi=a+b*(C+(1-C)*Z[:,3])
    Ei=hi-W-extra
    conf=hi>W
    cut=hi[conf]-W
    w0=hi[conf]-y
    s=np.minimum(cut/(xr*xr*w0),1.)
    pr=s/(1+np.sqrt(1-s))
    kin=np.maximum(2*a[conf]-hi[conf]-y,0.)
    rate=G[conf]*(1-C[conf])*(1+occupied(y,kin,d,W))*pr
    bind=Ei[conf]<0.
    delays=[]
    for h in TARGETS:
        # Exterior excursion <= full orbital period. Restrict sampling
        # to the energy band that can contribute at this time threshold.
        ep=.5*(2*np.pi*rt*extra*R0/(h*V0*GYR))**(2/3)
        c_min=np.maximum(C,(W+extra-ep-a)/b)
        c_max=np.minimum(1.,(W+extra-a)/b)
        sel=c_min<c_max
        if not np.any(sel):
            delays.append(0.)
            continue
        width=c_max[sel]-c_min[sel]
        e=a[sel]+b[sel]*(c_min[sel]+width*Z[sel,3])
        eps=np.maximum(W+extra-e,1e-15)
        v2=e-y
        s=np.clip((e-W)/(xr*xr*v2),0.,1.)
        pp=s/(1+np.sqrt(1-s))
        mu=1-pp*(1-Z[sel,4])
        l2=2*(xr*rt)**2*v2*(1-mu*mu)
        hold=duration(rt,extra,eps,l2)
        kin=np.maximum(2*a[sel]-e-y,0.)
        weight=G[sel]*width*(1+occupied(y,kin,d,W))*pp
        delays.append(weight[hold>h].sum()/len(G))
    return np.r_[n[-1]/np.sqrt(2),rate.sum()/len(G),rate[bind].sum()/len(G),delays]


def run(p):
    inputs=json.loads(p.read_text())['outputs']
    clock=V0**3*SIGMA_OVER_M*GYR/(4*np.pi*constants.G*R0**2)
    out=[]
    for row in inputs:
        F,W=row['peak_occupation'],row['escape_depth']
        d,h,n,N,sol,lr=solve(F,W,return_profile=True)
        extra=sol(lr)[1]
        rates=[]
        for nr,k,S in ((80,14,414),(112,17,414),(112,17,854)):
            x,w=roots_legendre(nr)
            z=-13+(lr+13)*(1+x)/2
            y=np.minimum(sol(z)[0],W-1e-13)
            sam=qmc.Sobol(5,scramble=True,seed=S).random_base2(k)
            a=np.array([local(s,d,W,np.exp(t-lr),np.exp(lr),extra,sam)
                        for s,t in zip(y,z,strict=True)])
            w3=w*np.exp(3*z)*(lr+13)/2
            assert abs((w3@a[:,0])/N-1)<8e-5
            rates.append(clock*(w3*a[:,0]**2)@a[:,1:]/(2*N*n))
        low,a,b=rates
        m=(a+b)/2
        assert abs(m[0]/row['actual_outward_accessible_rate']-1)<.012
        assert abs((m[1]/m[0])/row['negative_energy_among_outward_accessible_fraction']-1)<.016
        drift=np.max(abs(a/b-1))
        assert drift<.072,(W,rates)
        out.append(dict(escape_depth=W,peak_occupation=F,
                        outward_accessible_rate_per_Gyr=m[0],
                        outward_accessible_but_negative_energy_rate_per_Gyr=m[1],
                        fraction_of_negative_energy_outward_flux_exterior_longer_than_target=(m[2:]/m[1]).tolist(),
                        finest_scramble_maximum_relative_drift=float(drift),
                        coarse_to_fine_maximum_relative_difference=float(np.max(abs(low/m-1)))))
    return dict(scope="Ballistic exterior part only, fixed spherical vacuum exterior and no rescattering "
                      "in the time between ejection and exterior crossing. Return delay is conditional "
                      "on negative orbital energy and reaching the boundary. Not a formed debris halo.",
                exterior_duration_targets_Gyr=TARGETS,independent_radial_integral_defect=control(),outputs=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('crossing_ledger',type=Path)
    print(json.dumps(run(p.parse_args().crossing_ledger),indent=2,allow_nan=False))
