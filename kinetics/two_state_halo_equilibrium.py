"""Two-state near-degenerate Bose halo thermodynamics and initial conversion.

Both inertial/gravitational masses set equal (relative gap <1e-9), internal
excitation energy retained in collisions and chemical equilibrium. This is
a local box + isothermal Poisson comparison, NOT a new rotation curve fit.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.special import roots_legendre,zeta,gamma
from scipy.optimize import brentq,root
from scipy.integrate import solve_ivp
from soft_mixing_activation import scales,frac_B
from late_activation_overlap import critical_density
from resonant_mediator_inverse import M,C,SEC_GYR

ZS=float(zeta(1.5,1.))
Z5=float(zeta(2.5,1.))
coef={s:np.array([float(mp.zeta(s-j))*(-1)**j/math.factorial(j) for j in range(21)]) for s in [1.5,2.5]}
def pl(s,a):
    a=np.asarray(a,float);out=np.empty_like(a);b=a<.5
    out[b]=gamma(1-s)*a[b]**(s-1)+np.polynomial.polynomial.polyval(a[b],coef[s])
    x=a[~b];j=np.arange(1,81,dtype=float)
    out[~b]=np.sum(np.exp(-x[...,None]*j)/j**s,axis=-1)
    return out

def static_profiles(d,rmax=12.):
    # Original m,R,V0 fix rho(0) and T; choose the chemical potential
    # so the total degeneracy is the original one-state critical density.
    alpha=brentq(lambda a:float(pl(1.5,a)+pl(1.5,a+d)-ZS),0.,10.,xtol=1e-13)
    rr=np.geomspace(1e-5,rmax,160)
    ans=[]
    for a,two in [(0.,False),(alpha,True),(0.,True)]:
        n0=float(pl(1.5,a)+(pl(1.5,a+d) if two else 0))
        def rho(y):
            y=max(float(y),0.)
            return float(pl(1.5,a+y)+(pl(1.5,a+d+y) if two else 0))/n0
        eps=rr[0]
        def rhs(x,y):return [y[1]/x**2,x*x*rho(y[0])]
        sol=solve_ivp(rhs,(eps,rmax),[eps**2/6,eps**3/3],t_eval=rr,rtol=2e-9,atol=1e-12)
        assert sol.success
        dens=np.array([rho(y) for y in sol.y[0]])
        vel=np.sqrt(sol.y[1]/rr)
        ans.append(dict(two_state=two,alpha=a,nu_total=n0,
            central_heavy_fraction=float(pl(1.5,a+d)/n0) if two else 0.,
            density=dens.tolist(),circular_speed=vel.tolist(),
            central_slope_estimate=float((1-dens[0])/eps)))
    return dict(alpha_fixed_old_density=alpha,xi=rr.tolist(),profiles=ans,
        critical_density_ratio=float((ZS+pl(1.5,d))/ZS),
        critical_cusp_slope_ratio=float(ZS/(ZS+pl(1.5,d))),
        fixed_density_circular_speed_max_relative_change=float(max(abs(np.array(ans[1]['circular_speed'])/ans[0]['circular_speed']-1))))

def box_equilibrium(d):
    # Initial one-state ideal critical gas. Units T0, n0=ZS in common
    # phase-space normalization. Total energy = 3/2 Z5.
    E0=1.5*Z5
    def number(t,a):return t**1.5*(pl(1.5,a)+pl(1.5,a+d/t))
    def energy(t,a):return 1.5*t**2.5*(pl(2.5,a)+pl(2.5,a+d/t))+d*t**1.5*pl(1.5,a+d/t)
    tc=brentq(lambda t:float(energy(t,0.)-E0),.05,2.)
    if number(tc,0.)<ZS:
        a=0.;t=tc;cond=1-float(number(tc,0.))/ZS
    else:
        def alpha(t):return brentq(lambda a:float(number(t,a)-ZS),0.,30.,xtol=1e-13)
        tlo=brentq(lambda t:float(number(t,0.)-ZS),.05,2.)
        t=brentq(lambda t:float(energy(t,alpha(t))-E0),tlo*(1+1e-11),2.)
        a=alpha(t);cond=0.
    return dict(T_over_initial=t,alpha=a,condensate_number_fraction=cond,
        heavy_number_fraction=float(t**1.5*pl(1.5,a+d/t)/ZS),
        energy_relative_residual=float(energy(t,a)/E0-1),
        number_relative_residual=float(number(t,a)/ZS+cond-1))

def source(d,u,n=80,na=64,alpha=1e-10):
    # Initial critical L gas, H empty. Dimensionless speeds q=p/sqrt(2mT).
    # The relative-speed cross section convention is the measured LL sigma.
    q,w=roots_legendre(n);q=(q+1)*5;w=w*5
    w=w*q*q/np.expm1(q*q+alpha)
    N=w.sum();w/=N
    a,b=np.meshgrid(q,q,indexing='ij');wt=w[:,None]*w[None,:]
    # Angle integral using uniform relative cosine and threshold split.
    x,v=roots_legendre(na)
    totals=[]
    for dh,ratio,mixed in [(d,2*(1-u)/u,True),(2*d,((1-u)/u)**2,False)]:
        # energy q1^2+q2^2, relative k^2=|q1-q2|^2/4,
        # threshold k^2 >= dh/2.
        upper=np.minimum(1.,(a*a+b*b-2*dh)/(2*a*b))
        upper=np.maximum(upper,-1.)
        W=(upper+1)/2
        rate=np.zeros_like(a)
        for z,vv in zip(x,v):
            c=-1.+W*(z+1)
            P=np.sqrt(np.maximum(a*a+b*b+2*a*b*c,0.))
            k=np.sqrt(np.maximum((a*a+b*b-2*a*b*c)/4,0.))
            kf=np.sqrt(np.maximum(k*k-dh/2,0.))
            # vrel/vT=2k; final state flux ratio kf/k.
            enh=np.ones_like(P)
            if mixed:
                low=(P/2-kf)**2+alpha;high=(P/2+kf)**2+alpha
                span=2*P*kf
                # analytic average of f_L over final relative direction
                mid=(low+high)/2
                vsmall=span<1e-10*np.maximum(mid,1e-10)
                fmean=np.empty_like(P)
                fmean[vsmall]=1/np.expm1(mid[vsmall])
                fmean[~vsmall]=(np.log(-np.expm1(-high[~vsmall]))-np.log(-np.expm1(-low[~vsmall])))/span[~vsmall]
                enh+=np.maximum(fmean,0.)
            rate+=vv*W/2*2*kf*enh
        totals.append(float(ratio*np.sum(wt*rate)))
    # Pair events per L particle = .5*n_L*sigma_LL*vT * expectation.
    # LL->LH creates one H, LL->HH creates two.
    return dict(n=n,angle_order=na,alpha=alpha,normalized_radial_integral=float(N),
        LL_to_LH_event_kernel=totals[0],LL_to_HH_event_kernel=totals[1],
        heavy_initial_number_gain_over_n_sigma_vT=.5*totals[0]+totals[1])

def screen():
    p=scales(.1,z_on=5.,dneff=.1)
    u=p['u'];gap=p['mass_gap_eV']
    out=[]
    for V0 in [10.,20.,44.1,100.]:
        vT=np.sqrt(2)*V0*1e5;T=M*(vT/C)**2/2;d=gap/T
        t=box_equilibrium(d)
        c=[source(d,u,n,na) for n,na in [(48,48),(96,96),(144,144)]]
        # Use original critical L density and target LL sigma.
        dens=critical_density(vT)
        gain10=dens*p['sigma']*vT*10*SEC_GYR*c[-1]['heavy_initial_number_gain_over_n_sigma_vT']
        out.append(dict(V0_kms=V0,T_eV=T,gap_over_T=d,critical_density_cgs=dens,
             profiles=static_profiles(d),local_box_equilibrium=t,initial_conversion=c,
             initial_frozen_heavy_gain_10Gyr=gain10,
             scope='Fixed initial L critical bath. Frozen initial heavy source is not finite-time fraction. Equilibrium assumes chemical/kinetic relaxation, conserved number and energy, and ideal gas.'))
    return dict(params=p,rows=out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    r=screen()
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    for s in r['rows']:
        print(s['V0_kms'],s['gap_over_T'],s['local_box_equilibrium'],s['profiles']['alpha_fixed_old_density'],s['profiles']['fixed_density_circular_speed_max_relative_change'],s['initial_frozen_heavy_gain_10Gyr'],[x['heavy_initial_number_gain_over_n_sigma_vT'] for x in s['initial_conversion']],flush=True)
