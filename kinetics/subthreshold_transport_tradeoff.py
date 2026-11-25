"""Subthreshold virtual scalar mediator: transport, conversion, and loop tradeoff.

V includes (g/2) S B^2 and stabilizing lambda0 B^4/24 with
lambda0>=3g^2/mS^2. A/B softly mix at late times. The s-channel pole lies
below the light-light threshold; no instantaneous on-shell resonance.
In NR halo kinematics and for detuning <<m^2, retain the dominant s-channel
weight w_a(v)=[a/(a+(v/vT)^2)]^2, a=(4mL^2-mS^2)/(m^2 vT^2).
Compare a radial kurtosis Rayleigh rate, not total scattering count.
This is a conditional parameter screen, not a UV completion or fit.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,roots_laguerre,zeta
from scipy.integrate import quad,trapezoid
from scipy.interpolate import PchipInterpolator
from halo_bose_rayleigh import moments,quotient
from quartic_consistency import lambda_from_sigma
from direct_activation_consistency import rho_tot,radiation_T4,T0
from late_activation_overlap import H,critical_density
from soft_mixing_activation import frac_B
from resonant_mediator_inverse import M,C,SEC_GYR
from velocity_only_exposure import RHO0

def angle_integral(S,pmin,pmax,a):
    # P=total incoming momentum in units sqrt(2mT); vrel^2/vT^2=2S-P^2.
    # Exact primitive of a^2/(a+2S-P^2)^2.
    D=a+2*S;root=np.sqrt(D)
    def primitive(p):
        den=a+np.maximum(2*S-p*p,0.)
        x=np.minimum(p/root,np.nextafter(1.,0.))
        return .5*a*a*(p/(D*den)+np.arctanh(x)/(D*root))
    return np.maximum(primitive(pmax)-primitive(pmin),0.)

def rayleigh_grid(aa,alpha=1e-6,order=64):
    n,norm,A,B=moments(alpha)
    x,w=roots_legendre(order);u=(x+1)/2;w=w/2
    sl=np.log1p(80/alpha);S=alpha*np.expm1(sl*u);sw=w*(alpha+S)*sl
    out=np.zeros(len(aa));unweighted=0.
    for s,ds in zip(S,sw):
        al=np.log1p(.5*s/alpha)
        a=alpha/s*np.expm1(al*u);aw=w*(alpha/s+a)*al
        bl=np.log((alpha/s+.5)/(alpha/s+a))
        b=(alpha/s+a[:,None])*np.expm1(bl[:,None]*u)+a[:,None]
        bw=w*(alpha/s+b)*bl[:,None]
        f1=1/np.expm1(alpha+s*a[:,None]);f2=1/np.expm1(alpha+s*(1-a[:,None]))
        f3=1/np.expm1(alpha+s*b);f4=1/np.expm1(alpha+s*(1-b))
        diff=2*s*s*((b*b-b)-(a[:,None]**2-a[:,None]))
        base=ds*s*s*8*aw[:,None]*bw*f1*f2*(1+f3)*(1+f4)*diff**2
        weight=base.sum(axis=1)
        p1=np.sqrt(s*a);p2=np.sqrt(s*(1-a));lo=p2-p1;hi=p2+p1
        unweighted+=float(weight@(hi-lo))
        for j,d in enumerate(aa):
            out[j]+=weight@angle_integral(s,lo,hi,d)
    return out/(8*n*norm),unweighted/(8*n*norm)

def mb_integral(d,a,order=48):
    # Relative Maxwell speed w, threshold w^2>=2d. Incoming normalization
    # to MB particle number. Laguerre integration after exact threshold shift.
    y,wt=roots_laguerre(order)
    dd=np.asarray(d)[...,None]
    return 2*np.sqrt(2/np.pi)*np.exp(-np.minimum(np.asarray(d),745))*np.sum(
        wt*np.sqrt(y*(y+dd))*(a/(a+2*(y+dd)))**2,axis=-1)

def bose_source(d,u,a,n=96,na=96,alpha=1e-6):
    q,w=roots_legendre(n);q=(q+1)*5;w=w*5
    w=w*q*q/np.expm1(q*q+alpha);w/=w.sum()
    p,b=np.meshgrid(q,q,indexing='ij');wt=w[:,None]*w[None,:]
    x,v=roots_legendre(na);totals=[]
    for dh,ratio,mixed in [(d,2*(1-u)/u,True),(2*d,((1-u)/u)**2,False)]:
        upper=np.clip((p*p+b*b-2*dh)/(2*p*b),-1.,1.);W=(upper+1)/2
        rate=np.zeros_like(p)
        for z,vv in zip(x,v):
            c=-1.+W*(z+1)
            P=np.sqrt(np.maximum(p*p+b*b+2*p*b*c,0.))
            k2=np.maximum((p*p+b*b-2*p*b*c)/4,0.)
            kf=np.sqrt(np.maximum(k2-dh/2,0.))
            enh=np.ones_like(P)
            if mixed:
                low=(P/2-kf)**2+alpha;high=(P/2+kf)**2+alpha
                span=2*P*kf;mid=(low+high)/2
                small=span<1e-10*np.maximum(mid,1e-10);fm=np.empty_like(P)
                fm[small]=1/np.expm1(mid[small])
                fm[~small]=(np.log(-np.expm1(-high[~small]))-np.log(-np.expm1(-low[~small])))/span[~small]
                enh+=np.maximum(fm,0.)
            rate+=vv*W/2*2*kf*enh*(a/(a+4*k2))**2
        totals.append(float(ratio*np.sum(wt*rate)))
    return .5*totals[0]+totals[1]

def mass_tradeoff(a,Q,V,z,dneff,alpha=1e-6,efold=1.,age=10.,threshold_factor=1.):
    # Optimistically assume full activation at fixed halo conditions for
    # all 10 Gyr. For actual activation, required coupling is larger.
    vT=np.sqrt(2)*V*1e5;T=M*(vT/C)**2/2
    rho=critical_density(vT);scale=2*M*T # eV^2
    det=a*scale;ms2=4*M*M-det
    if ms2<=0:raise ValueError('mediator mass squared not positive')
    sigma0=efold/(rho*vT*age*SEC_GYR*Q)
    lam0=float(lambda_from_sigma(sigma0))
    # Dominant s-channel normalization. Stabilizing contact is included
    # in the illustrative loop floor but not used to artificially increase
    # the halo relaxation rate.
    U=np.geomspace(1e-9,.49,1800)
    g2=lam0*det/U**2;quartic=3*g2/ms2
    Dfloor=threshold_factor*(g2/(16*np.pi**2)+quartic*M*M/(32*np.pi**2))
    cap=min(np.pi**2/9*radiation_T4(z,dneff),.04*rho_tot(z))
    delta_max=4*np.pi*np.sqrt(cap)
    Dmax=delta_max*(1-2*U)/np.sqrt(U*(1-U))
    dmax=delta_max/(scale*np.sqrt(U*(1-U)))
    # Lower bound on Bose initial conversion: replace incoming
    # f>=exp(-alpha-q^2), final-state stimulation>=1. Normalization to
    # BE density gives ratio exp(-alpha)^2/Li_{3/2}(exp(-alpha))^2.
    from two_state_halo_equilibrium import pl
    ratio=(np.exp(-alpha)/float(pl(1.5,alpha)))**2
    r=(1-U)/U
    mb=(r*mb_integral(dmax,a)+r*r*mb_integral(2*dmax,a))*ratio/Q
    # Keep separate: splitting threshold only; perturbativity and
    # illustrative mediator pole stability (1-loop threshold <= detuning).
    stable=g2/(32*np.pi**2)<=det
    masks=[Dfloor<=Dmax,(Dfloor<=Dmax)&stable,(Dfloor<=Dmax)&stable&(quartic<1.)&(g2/ms2<1.)]
    result=[]
    for label,ok in zip(['splitting_only','also_pole_threshold','also_perturbative'],masks):
        if not np.any(ok):result.append(dict(screen=label,feasible=False));continue
        i=np.flatnonzero(ok)[np.argmin(mb[ok])]
        result.append(dict(screen=label,feasible=True,u=float(U[i]),gap_over_T=float(dmax[i]),
           lower_Bose_conversion_per_local_efold=float(mb[i]),D_floor_over_cap=float(Dfloor[i]/Dmax[i]),
           g_eV=float(np.sqrt(g2[i])),lambda_stabilizer=float(quartic[i]),
           mediator_loop_over_detuning=float(g2[i]/(32*np.pi**2*det)),
           D_eV2=float(Dmax[i]),delta_eV2=float(delta_max),
           vacuum_threshold_over_cap=1.))
    return dict(a=a,Q=Q,V0_kms=V,z_on=z,dneff=dneff,efold=efold,sigma0_needed=sigma0,
       kinetic_scale_eV2=scale,detuning_eV2=det,results=result,
       scope='NR s-channel and fixed 10 Gyr halo. Illustrative loop scales, not invariant naturalness bounds. MB replacement is a lower bound on initial Bose conversion within this collision model.')

def cosmic_proxy(row,case,aa,Qgrid,refQ):
    p=case['params'];h=case['result']['history']
    r=np.array([v['z']+1 for v in h]);x=np.array([v['x'] for v in h]);_,ix=np.unique(r,return_index=True);ix=np.sort(ix);r=r[ix];x=x[ix]
    late=r<=101;r,x=r[late],x[late];t=-np.log(r)
    Tr=T0*(3.909/row['g_bath'])**(1/3)/row['final']['theta']
    v0=row['final']['peak_bin_q']*Tr/M*C
    vT=np.sqrt(2)*44.1e5;scale=M*M*(vT/C)**2
    d=p['delta'];D=p['Delta'];u=p['u']
    shift=lambda x:-2*d*d*x*x/(np.hypot(D,2*d*x)+D)
    change=4*(shift(x)-shift(1.))/scale
    Qfun=PchipInterpolator(np.log(aa),np.log(Qgrid),extrapolate=False)
    out=[]
    rho=critical_density(vT)
    for a in aa:
        aval=a+change
        if min(aval)<=0.:
            out.append(dict(a_at_full=a,below_threshold_all_late_times=False,minimum_a=float(min(aval))))
            continue
        mix=frac_B(x,d,D)
        shape=(mix/u)**4*(a/aval)**2
        Q=np.exp(Qfun(np.log(np.clip(aval,aa[0],aa[-1]))))
        halo=rho*vT*trapezoid(shape*Q/H(r-1),t)
        sig=1/halo
        tau=sig*trapezoid(RHO0*v0*r**4*shape/H(r-1),t)
        out.append(dict(a_at_full=a,below_threshold_all_late_times=True,minimum_a=float(min(aval)),
           sigma0_at_full_for_one_efold=float(sig),late_tau=float(tau),
           frozen_relative_Ngt2_proxy=float(-5.756e7*tau),
           maximum_cold_relative_speed_term_over_detuning=float(max((12*v0*r/vT)**2/aval)),
           minimum_virtual_lifetime_s=float(6.582119569e-16*4*M/(scale*max(aval))),
           maximum_virtual_lifetime_s=float(6.582119569e-16*4*M/(scale*min(aval)))))
    return dict(params=p,scheme=case['result']['loop_scheme'],rows=out,
       scope='Prescribed original soft-mixing roll, rescale scattering normalization for one local Rayleigh e-fold. Approximate cosmic contact kernel only if cold speed term small; frozen moment proxy not actual nonlinear change. Halo fixed density all times.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--cases',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0]
    aa=np.geomspace(1e-5,1e3,81);Q1,Q0=rayleigh_grid(aa,order=64);Q2,Q00=rayleigh_grid(aa,order=96)
    err=float(max(abs(Q2/Q1-1)))
    rows=[]
    for V in [10.,44.1,100.]:
        for z in [3.,5.,8.,10.,20.]:
            for dneff in [.1,.3]:
                for av,Q in zip(aa,Q2):
                    rows.append(mass_tradeoff(float(av),float(Q),V,z,dneff))
    best=[]
    for V in [10.,44.1,100.]:
        for z in [3.,5.,8.,10.,20.]:
            for dneff in [.1,.3]:
                for lab in ['splitting_only','also_pole_threshold','also_perturbative']:
                    arows=[(r,x) for r in rows if r['V0_kms']==V and r['z_on']==z and r['dneff']==dneff for x in r['results'] if x['screen']==lab and x['feasible']]
                    if not arows:continue
                    r,x=min(arows,key=lambda z:z[1]['lower_Bose_conversion_per_local_efold'])
                    best.append(dict(a=r['a'],Q=r['Q'],V0_kms=V,z_on=z,dneff=dneff,sigma0_needed=r['sigma0_needed'],**x))
    cases=json.loads(Path(a.cases).read_text())
    cosmic=[cosmic_proxy(row,c,aa,Q2,Q00) for c in cases if c['params']['sigma'] in [.1,.3,1.25]]
    # Independent exact Bose evaluation of the most optimistic F583 points.
    exact=[]
    for b in best:
        if b['V0_kms']!=44.1 or b['z_on'] not in [5.,8.,20.] or b['dneff']!=.1 or b['screen']!='also_perturbative':continue
        val=[bose_source(b['gap_over_T'],b['u'],b['a'],n,n) for n in [48,96,144]]
        exact.append(dict(**b,bose_conversion_per_efold=[v/b['Q'] for v in val]))
    # Constant-cross-section and direct angular checks.
    checks=[]
    for s,l,h,av in [(1.,.3,1.2,.001),(2.,.01,.02,1e-5),(3.,.1,2.,10.),(2.,.1,1.9,1000.)]:
        x=quad(lambda P:(av/(av+2*s-P*P))**2,l,h,epsabs=1e-20,epsrel=1e-10)[0]
        y=float(angle_integral(s,l,h,av));checks.append(dict(S=s,a=av,relative_error=y/x-1))
    out=dict(aa=aa.tolist(),Q64=Q1.tolist(),Q96=Q2.tolist(),contact64=Q0,contact96=Q00,
      contact_existing=quotient(np.exp(-1e-6),96)['rayleigh_over_n_sigma_vT'],max_grid_relative_difference=err,
      angle_checks=checks,rows=rows,best=best,exact=exact,cosmic=cosmic,
      scope='Subthreshold s-channel optimistic screen. Fixed halo Rayleigh mode is not global equilibration. Loop floor/cap model dependent. No cosmological nonlinear evolution, UV matching, or mediator population.')
    Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('Q',Q00,'max err',err,'exact',checks)
    for b in exact:print('BEST',b,flush=True)
    for c in cosmic:
        v=[v for v in c['rows'] if v['below_threshold_all_late_times'] and abs(v['frozen_relative_Ngt2_proxy'])<.01]
        print('COSMIC',c['params']['sigma'],c['scheme'],'min a under 1% proxy',min([r['a_at_full'] for r in v],default=None),flush=True)

if __name__=='__main__':main()
