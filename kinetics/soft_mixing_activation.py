"""Conditional softly mixed two-state late-activation model.

V = m^2 A^2/2 + (m^2+Delta) B^2/2 + kappa chi A B
  + lambda_B B^4/24 + lambda_chi (chi^2-v^2)^2/4
  + g^2 chi^2 sum eta_i^2/2.
At chi=0 the lighter state A has no tree B quartic. Trace of mass matrix
is chi independent. This is a low-energy parameter screen, not a UV
completion or demonstrated cosmological/halo solution.
"""
import argparse,json,time,math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,trapezoid,quad
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre,zeta
from quartic_consistency import lambda_from_sigma,shell_nodes
from late_activation_overlap import H,H0,Om,Or,critical_density
from direct_activation_consistency import radiation_T4,rho_tot,CGS_RHO_TO_EV4,T0,MPL_REDUCED
from resonant_mediator_inverse import HBAR,C,M
from velocity_only_exposure import RHO0,SEC_GYR,halo_pdf

UOPT=(3-np.sqrt(3))/4

def frac_B(x,delta,Delta):
    # Rationalized form accurate for a tiny mixing.
    r=2*delta*x/Delta;s=np.sqrt(1+r*r)
    return r*r/(2*s*(s+1))

def scales(sigma=.1,u=UOPT,z_on=5.,dneff=.1,split=1.,X=None,muH=30.,Neta=1):
    lam=float(lambda_from_sigma(sigma));lamB=lam/(u*u)
    Delta=split*lamB*M*M/(32*np.pi**2)
    delta=Delta*np.sqrt(u*(1-u))/(1-2*u)
    # Chi-dependent trace-free 1-loop scale; exact MS finite piece can
    # vanish at a chosen matching scale. This is a log-threshold diagnostic.
    Xmin=delta*delta/(16*np.pi**2)
    Xstable=6*delta**4/(M**4*lamB)
    Td4=radiation_T4(z_on,dneff)
    Xmax=min(np.pi**2/9*Td4,.04*rho_tot(z_on))
    if X is None:X=np.sqrt(Xmin*Xmax)
    mu=muH*H(z_on)*HBAR;v=np.sqrt(X)/mu
    Td=(Td4/Neta)**.25
    m_eta=np.sqrt(12*X/(Neta*Td*Td))
    kap=delta/v;g=m_eta/v;lamchi=X/v**4
    dm2=np.sqrt(Delta*Delta+4*delta*delta)
    dm=dm2/(np.sqrt(M*M+(Delta+dm2)/2)+np.sqrt(M*M+(Delta-dm2)/2))
    shift2=-2*delta*delta/(np.sqrt(Delta*Delta+4*delta*delta)+Delta)
    mL=np.sqrt(M*M+shift2)
    # Derivative of m_L wrt chi (dimensionless Yukawa-like scalar charge).
    charge=-delta*delta/(mL*v*dm2)
    fifth=2*(MPL_REDUCED*charge/mL)**2
    wall=(2*np.sqrt(2)/3)*mu*v*v
    return dict(sigma=sigma,u=u,z_on=z_on,dneff=dneff,split=split,mu_over_H_on=muH,Neta=Neta,
        lambda_late=lam,lambda_B=lamB,Delta=Delta,delta=delta,X=X,X_loop_floor=Xmin,X_stability_floor=Xstable,X_thermal_or_vacuum_cap=Xmax,
        overlap_ratio=Xmax/max(Xmin,Xstable),passes_illustrative_loop_screen=bool(max(Xmin,Xstable)<=X<=Xmax),
        mu_eV=mu,vev_eV=v,vev_over_Mpl=v/MPL_REDUCED,kappa_eV=kap,g_eta=g,lambda_chi=lamchi,
        Td_on_eV=Td,m_eta_at_vev_eV=m_eta,eta_mass_over_Td_on=m_eta/Td,
        mass_gap_eV=dm,light_mass_shift_relative=shift2/(M*(mL+M)),fifth_force_over_gravity=fifth,
        vacuum_energy_over_total=X/4/rho_tot(z_on),wall_scaling_fraction_on=wall*HBAR*H(z_on)/rho_tot(z_on),
        eta_relic_density_fraction_NR=Neta*zeta(3,1)/np.pi**2*(Td/(1+z_on))**3*m_eta/(RHO0*CGS_RHO_TO_EV4))

def bath_integrals(refined=True):
    if refined:
        u,w=roots_legendre(32);e=np.array([0.,1e-5,1e-4,1e-3,.01,.1,1.,5.,15.,50.])
        y=((e[:-1,None]+e[1:,None])/2+np.diff(e)[:,None]*u/2).ravel()
        w=(np.diff(e)[:,None]*w/2).ravel()
    else:
        u,w=roots_legendre(96);y=(u+1)*25;w=w*25
    f=1/np.expm1(y)
    j=np.geomspace(1e-6,1e5,600)
    variance=np.sum(w[None,:]*y[None,:]**2*f[None,:]/np.sqrt(y[None,:]**2+j[:,None]**2),axis=1)/(2*np.pi**2)
    interp=PchipInterpolator(np.log(j),np.log(variance))
    def V(x):
        x=np.asarray(x)
        return np.where(x<=j[0],1/12,np.where(x>=j[-1],zeta(3,1)/(np.pi**2*np.maximum(x,1e-300)),np.exp(interp(np.log(np.clip(x,j[0],j[-1]))))))
    return V

def variance_table(row):
    e=np.array(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    q,w=shell_nodes(e,f,12);w/=w.sum()
    theta=row['final']['theta'];Tref=T0*(3.909/row['g_bath'])**(1/3)/theta
    r=np.geomspace(1.,1e16,1200)
    inv=np.sum(w[None,:]/np.hypot(M,r[:,None]*Tref*q[None,:]),axis=1)
    N=RHO0*CGS_RHO_TO_EV4/M
    v=N*r**3*inv
    interp=PchipInterpolator(np.log(r),np.log(v))
    return lambda r:np.exp(interp(np.log(np.clip(r,1.,1e16)))),Tref

def halo_channels(p,x=1.,vT=np.sqrt(2)*44.1e5):
    u=frac_B(x,p['delta'],p['Delta'])
    dm2=np.hypot(p['Delta'],2*p['delta']*x);dm=dm2/(2*M)
    t1=C*np.sqrt(4*dm/M);t2=C*np.sqrt(8*dm/M)
    # Low-energy near-degenerate equal-mass approximation, distinguishable
    # LH final state factor two. Initial LL only, not a solved two-species gas.
    def avg(t):
        return quad(lambda z:z*halo_pdf(z,vT)*np.sqrt(max(0.,1-(t/z)**2)),
                    t,max(t,12*vT),epsabs=1e-5,epsrel=1e-8)[0]
    base=2*np.sqrt(2/np.pi)*vT
    r1=2*(1-u)/max(u,1e-300)*avg(t1)/base
    r2=((1-u)/max(u,1e-300))**2*avg(t2)/base
    T=M*(vT/C)**2/2
    # Polylog Li_3/2(e^-dm/T) via convergent integral (normalized).
    d=dm/T
    li=quad(lambda t:2*t*t/np.expm1(t*t+d),0,10,epsrel=1e-9,epsabs=1e-11)[0]/(np.sqrt(np.pi)/2)
    gcrit=1+li/zeta(1.5,1)
    return dict(x=x,u=float(u),mass_gap_eV=float(dm),vT_kms=vT/1e5,
        threshold_LH_kms=t1/1e5,threshold_HH_kms=t2/1e5,
        rate_LH_over_LL=float(r1),rate_HH_over_LL=float(r2),
        critical_density_ratio_two_species=float(gcrit),illustrative_mass_at_fixed_critical_density_eV=M/gcrit**.25)

def evolve(row,p,seed=1e-8,include_medium=True,rtol=2e-7,loop_scheme='none',refined_bath=True):
    start=time.monotonic();V,Tref=variance_table(row);etaV=bath_integrals(refined_bath)
    mu=p['mu_eV'];X=p['X'];Delta=p['Delta'];delta=p['delta'];r_on=1+p['z_on'];u=p['u']
    eta0=p['m_eta_at_vev_eV'];td0=p['Td_on_eV']/r_on
    # Renormalized eta Coleman-Weinberg contribution. Two explicit
    # matching prescriptions: fix origin curvature and vev, or vev
    # and its curvature. A constant vacuum offset is immaterial to roll.
    acw=p['Neta']*eta0**4/(64*np.pi**2)
    def loop_term(x):
        x2=x*x;logx=np.log(max(x2,1e-300))
        if loop_scheme=='origin-vev':return 4*acw/X*x2*logx
        if loop_scheme=='vev-curvature':return 4*acw/X*(1+x2*(logx-1))
        if loop_scheme=='none':return 0.
        raise ValueError(loop_scheme)
    def loop_energy(x):
        x2=x*x;logx=np.log(max(x2,1e-300))
        if loop_scheme=='origin-vev':return acw*(x2*x2*(logx-.5)+.5)
        if loop_scheme=='vev-curvature':return acw*(x2*x2*(logx-1.5)+2*x2-.5)
        return 0.

    def fun(t,y):
        r=np.exp(-t);x,dx=y;hh=H(r-1)*HBAR
        dH=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        therm=12*(r/r_on)**2*etaV(eta0*abs(x)/(td0*r))
        med=delta*delta*V(r)/(X*np.hypot(Delta,2*delta*x)) if include_medium else 0.
        return [dx,-(3+dH)*dx-(mu/hh)**2*x*(x*x-1+therm-med+loop_term(x))]
    # Earlier frozen PSD extrapolation gives an adiabatic ground-state
    # medium susceptibility upper estimate, not a dynamical flavor calculation.
    t=np.linspace(-np.log(1e16),0.,2400)
    sol=solve_ivp(fun,(t[0],0.),[seed,0.],t_eval=t,method='DOP853',rtol=rtol,atol=rtol*max(seed,1e-8)*.01,max_step=.04)
    assert sol.success,sol.message
    r=np.exp(-t);x=sol.y[0];dx=sol.y[1];u_t=frac_B(x,delta,Delta)
    sig=p['sigma']*(u_t/u)**4
    # A single low-density tree 2->2 total cross-section upper count ignoring
    # kinematic threshold and stimulation: sum LL,LH,HH gives lambda_B^2*u^2.
    # Gamma <= amplitude^2*n*<1/E>^2/(64pi). This does not bound Bose evolution.
    n=RHO0*CGS_RHO_TO_EV4/M*r**3;var=V(r)
    dilute=p['lambda_B']**2*u_t**2*var**2/(64*np.pi*n*HBAR*H(r-1))
    # The initial heavy-state collision source HH ignores stimulation.
    HH=dilute*(1-u_t)**2
    # Use z>=1e3 as early window; integrating over all higher redshifts on
    # frozen PSD is intentionally conservative for total dilute count.
    early=(r>=1001)
    exp_early=float(trapezoid(HH[early],t[early]))
    mix_early=float(max(u_t[early]));sig_early=float(max(sig[early]))
    # Time varying contact diagnostic after z=100 (full NR distribution safe).
    late=r<=101;v0=row['final']['peak_bin_q']*Tref/M*C
    tau=float(trapezoid(RHO0*v0*r[late]**4*sig[late]/H(r[late]-1),t[late]))
    halo_t=float(trapezoid(sig/H(r-1),t))
    halo_late=float(trapezoid(sig[late]/H(r[late]-1),t[late]))
    rhoh=critical_density(np.sqrt(2)*44.1e5);vT=np.sqrt(2)*44.1e5
    # One previous trial radial Rayleigh decay, using only LL elastic.
    efold=.791935071*rhoh*vT*halo_late
    static=delta*delta*var/(Delta*X)
    thermo=(r/r_on)**2
    instability=np.maximum(1+static-thermo-loop_term(0.),0.)
    max_early_growth_upper=float(max((mu/(HBAR*H(r-1))*np.sqrt(instability))[early]))
    threshold={}
    for v in [.01,.1,.5,1.]:
        k=np.flatnonzero(sig>=v*p['sigma'])
        threshold[str(v)]=float(r[k[0]]-1) if len(k) else None
    # Today free-streaming bath contribution and chi homogeneous energy,
    # excluding the arbitrary cosmological constant offset.
    rhot=rho_tot(0.)
    field_energy=(.5*p['vev_eV']**2*(HBAR*H(0.)*dx[-1])**2+X/4*(x[-1]**2-1)**2+loop_energy(x[-1]))/rhot
    return dict(seed=seed,include_medium=include_medium,loop_scheme=loop_scheme,refined_bath=refined_bath,eta_CW_over_X=acw/X,zero_temperature_origin_mass_squared_over_mu2=-1+loop_term(0.),renormalized_false_vacuum_energy_over_total=(X/4+loop_energy(0.))/rho_tot(p['z_on']),rtol=rtol,wall_seconds=time.monotonic()-start,
        today_x=float(x[-1]),today_sigma_LL=float(sig[-1]),first_sigma_fraction_crossings=threshold,
        maximum_early_mixing_probability=mix_early,maximum_early_sigma_LL=sig_early,
        early_dilute_HH_count_upper=exp_early,maximum_early_adiabatic_tachyon_rate_over_H=max_early_growth_upper,
        integrated_sigma_seconds=halo_t,late_z100_frozen_NR_tau=tau,late_z100_frozen_N_gt2_relative_change_proxy=tau*(-5.756e7),
        conditional_local_LL_rayleigh_efolds=efold,field_energy_fraction_today=float(field_energy),
        halo=halo_channels(p,x[-1]),history=[dict(z=float(r[i]-1),x=float(x[i]),dx_dln_a=float(dx[i]),
                mixing_B=float(u_t[i]),sigma_LL=float(sig[i]),static_medium_term=float(static[i]),thermal_at_origin=float(thermo[i]),
                dilute_HH_upper_per_H=float(HH[i])) for i in range(0,len(t),8)]+[dict(z=0.,x=float(x[-1]),dx_dln_a=float(dx[-1]),mixing_B=float(u_t[-1]),sigma_LL=float(sig[-1]))],
        scope='Conditional homogeneous scalar roll, collisionless eta distribution, instantaneous light-eigenstate medium (upper adiabatic negative susceptibility), frozen primordial A PSD; NO quantum flavor conversion, full halo mixture, perturbations or UV protection.')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--input',required=True);a.add_argument('--output',required=True)
    a.add_argument('--sigma',type=float,default=.1);a.add_argument('--z-on',type=float,default=5.);a.add_argument('--dneff',type=float,default=.1)
    a.add_argument('--split',type=float,default=1.);a.add_argument('--X',type=float)
    a.add_argument('--muH',nargs='+',type=float,default=[10.,30.,100.]);a.add_argument('--seeds',nargs='+',type=float,default=[1e-8,1e-4])
    a.add_argument('--rtol',type=float,default=2e-7);args=a.parse_args()
    row=json.loads(Path(args.input).read_text())[0];res=[]
    for muH in args.muH:
        p=scales(args.sigma,UOPT,args.z_on,args.dneff,args.split,args.X,muH)
        for seed in args.seeds:
            r=evolve(row,p,seed,rtol=args.rtol);res.append(dict(params=p,result=r))
            Path(args.output).write_text(json.dumps(res,indent=2,allow_nan=False)+'\n')
            print('PARAMS',json.dumps(p),'RESULT',json.dumps({k:v for k,v in r.items() if k!='history'}),flush=True)
