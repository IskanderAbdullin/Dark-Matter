"""Time-dependent trigger scan for the virtual subthreshold mediator.

Use the existing homogeneous scalar/eta model, but do NOT identify
lambda_B with the low-energy resonant amplitude. Instead the stabilizing
quartic is 3 g^2/mS^2, and the small mass-splitting floor includes the
cubic one-loop threshold. For small mixing u, x=chi/v gives sigma_LL~x^8,
LL->LH~x^6, LL->HH~x^4. Use this only to shortlist, then exact roll,
detuning and mixing weights on the chosen candidates.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.special import zeta
from soft_mixing_activation import variance_table,bath_integrals,frac_B
from late_activation_overlap import H,Om,Or,critical_density
from direct_activation_consistency import radiation_T4,rho_tot,T0,MPL_REDUCED
from resonant_mediator_inverse import HBAR,M,C,SEC_GYR
from velocity_only_exposure import RHO0
from quartic_consistency import lambda_from_sigma
from subthreshold_transport_tradeoff import mb_integral,bose_source,rayleigh_grid
from two_state_halo_equilibrium import pl

def params(z,dneff,Xfrac,margin,muH,u=1e-5,Neta=1):
    T4=radiation_T4(z,dneff);cap=min(np.pi**2/9*T4,.04*rho_tot(z))
    X=Xfrac*cap;delta=4*np.pi*np.sqrt(X/margin)
    D=delta*(1-2*u)/np.sqrt(u*(1-u));mu=muH*H(z)*HBAR
    v=np.sqrt(X)/mu;Td=(T4/Neta)**.25;eta=np.sqrt(12*X/(Neta*Td*Td))
    return dict(z_on=z,dneff=dneff,Xfrac=Xfrac,loop_margin=margin,mu_over_H_on=muH,u=u,Neta=Neta,
       X=X,Xcap=cap,delta=delta,Delta=D,mu_eV=mu,vev_eV=v,Td_on_eV=Td,m_eta_at_vev_eV=eta,
       lambda_chi=X/v**4,eta_CW_over_X=Neta*eta**4/(64*np.pi**2*X),
       eta_mass_over_T=eta/Td,vev_over_Mpl=v/MPL_REDUCED,seed=1e-8)

def roll(row,p,scheme,seed=1e-8,rtol=4e-7,n=12000):
    V,Tr=variance_table(row);eta=bath_integrals()
    D=p['Delta'];delta=p['delta'];mu=p['mu_eV'];X=p['X'];r_on=p['z_on']+1
    m_eta=p['m_eta_at_vev_eV'];td0=p['Td_on_eV']/r_on;cw=p['eta_CW_over_X']
    def fun(t,y):
        r=np.exp(-t);x,dx=y;hh=HBAR*H(r-1)
        dh=-.5*(3*Om*r**3+4*Or*r**4)/(Om*r**3+Or*r**4+1-Om-Or)
        therm=12*(r/r_on)**2*eta(m_eta*abs(x)/(td0*r))
        med=delta**2*V(r)/(X*np.hypot(D,2*delta*x))
        xx=x*x;lx=np.log(max(xx,1e-300))
        cwterm=4*cw*(xx*lx if scheme=='origin-vev' else 1+xx*(lx-1))
        return [dx,-(3+dh)*dx-(mu/hh)**2*x*(xx-1+therm-med+cwterm)]
    # Start sufficiently before radiation-matter equality, check initial
    # frozen history against alternate start when refining selected cases.
    t0=-np.log(1e12)
    ts=np.r_[np.linspace(t0,-np.log(101.),600),np.linspace(-np.log(101.),0.,n)[1:]]
    sol=solve_ivp(fun,(t0,0.),[seed,0.],t_eval=ts,rtol=rtol,atol=rtol*1e-10,method='DOP853',max_step=.045)
    assert sol.success,sol.message
    r=np.exp(-ts);x=sol.y[0];dx=sol.y[1]
    v0=row['final']['peak_bin_q']*Tr/M*C
    return dict(t=ts,r=r,x=x,dx=dx,nfev=sol.nfev,v0=v0,
       max_early_abs_field=float(max(abs(x[r>=1001]))),today_field=float(x[-1]),
       z_at_half_x8=float(r[np.flatnonzero(abs(x)**8>=.5)[0]]-1) if np.any(abs(x)**8>=.5) else None,
       maximum_abs_x=float(max(abs(x))))

def moments_roll(rr,formation_z=2.):
    t,r,x=rr['t'],rr['r'],rr['x'];late=r<=101;halo=r<=formation_z+1
    clock=1/H(r-1)
    A={j:float(trapezoid(abs(x[halo])**(2*j)*clock[halo],t[halo])) for j in [2,3,4]}
    cosmic=float(trapezoid(r[late]**4*abs(x[late])**8*clock[late],t[late]))
    return A,cosmic

def shortlist(p,rr,aa,QQ,formation_z=2.,target_dist=.01,target_conv=.1):
    vT=np.sqrt(2)*44.1e5;rho=critical_density(vT);T=M*(vT/C)**2/2;scale=2*M*T
    A,cosmic=moments_roll(rr,formation_z)
    if A[4]<1e10:return dict(feasible=False,reason='no appreciable activation',A4_Gyr=A[4]/SEC_GYR)
    U=np.geomspace(1e-8,.49,1600)
    D=p['delta']*(1-2*U)/np.sqrt(U*(1-U));d=p['delta']/(scale*np.sqrt(U*(1-U)))
    ratio=(np.exp(-1e-6)/float(pl(1.5,1e-6)))**2
    rU=(1-U)/U
    results=[]
    for a,Q in zip(aa,QQ):
        sig=1/(rho*vT*Q*A[4]);lam=float(lambda_from_sigma(sig));det=a*scale
        g2=lam*det/U**2;ms2=4*M*M-det;quartic=3*g2/ms2
        floor=g2/(16*np.pi**2)+quartic*M*M/(32*np.pi**2)
        pscale=g2/(32*np.pi**2*det)
        cosmic_frac=5.756e7*sig*RHO0*rr['v0']*cosmic
        mb=ratio*(rU*mb_integral(d,a)*A[3]/A[4]+rU*rU*mb_integral(2*d,a)*A[2]/A[4])/Q
        # Keep all screened rows, including cosmic failures. Choose exact
        # candidates from those passing approximate cosmic proxy.
        for lab,mask in [('loop_thresholds',(floor<=D)&(pscale<=1)),
                         ('factor2_margin',(2*floor<=D)&(2*pscale<=1)),
                         ('factor10_margin',(10*floor<=D)&(10*pscale<=1))]:
            ok=mask&(quartic<1.)&(g2/ms2<1.)
            if not np.any(ok):continue
            i=np.flatnonzero(ok)[np.argmin(mb[ok])]
            results.append(dict(screen=lab,a=float(a),Q=float(Q),u=float(U[i]),D=float(D[i]),gap_over_T=float(d[i]),
              sigma0_for_one_efold=float(sig),g_eV=float(np.sqrt(g2[i])),lambda0=float(quartic[i]),
              split_floor_ratio=float(floor[i]/D[i]),pole_loop_ratio=float(pscale[i]),
              lower_conversion_per_efold=float(mb[i]),cosmic_frozen_proxy=float(cosmic_frac),
              proxy_pass=bool(cosmic_frac<=target_dist)))
    best={}
    for lab in ['loop_thresholds','factor2_margin','factor10_margin']:
        v=[r for r in results if r['screen']==lab and r['proxy_pass']]
        if v:best[lab]=min(v,key=lambda r:r['lower_conversion_per_efold'])
    return dict(feasible=bool(best),formation_z=formation_z,A2_Gyr=A[2]/SEC_GYR,A3_Gyr=A[3]/SEC_GYR,A4_Gyr=A[4]/SEC_GYR,
       cosmic_contact_factor=float(5.756e7*RHO0*rr['v0']*cosmic/(rho*vT*A[4])),
       best=best,rows=results)

def scan(row,aa,QQ,output,group='all'):
    # Focus on regimes where static threshold-only screen was promising.
    cases=[
      (8,.1,.7,1.1,20),(8,.1,.7,1.1,30),(8,.1,.9,1.1,20),
      (8,.2,.7,1.1,20),(8,.2,.7,1.1,30),(8,.3,.7,1.1,20),
      (10,.1,.7,1.1,20),(10,.1,.7,1.1,30),(10,.1,.5,2.,30),
      (10,.2,.5,1.1,20),(10,.2,.7,2.,30),(12,.1,.7,1.1,20),
      (12,.1,.5,2.,30),(12,.2,.5,2.,20),(15,.1,.7,1.1,15),
      (15,.1,.7,1.1,30),(20,.1,.7,1.1,10),(20,.1,.7,1.1,20),
      (20,.2,.5,2.,10),(20,.2,.5,2.,20)]
    out=[]
    for j,(z,dn,xf,margin,muH) in enumerate(cases):
        if group=='even' and j%2:continue
        if group=='odd' and not j%2:continue
        p=params(z,dn,xf,margin,muH)
        for scheme in ['origin-vev','vev-curvature']:
            start=time.monotonic();r=roll(row,p,scheme)
            q=shortlist(p,r,aa,QQ,2.)
            a=dict(params=p,scheme=scheme,roll={k:v for k,v in r.items() if k not in ['t','r','x','dx']},
                   screen=q,wall_seconds=time.monotonic()-start)
            out.append(a);Path(output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print('SCAN',z,dn,xf,margin,muH,scheme,r['z_at_half_x8'],
                {k:(v['a'],v['u'],v['gap_over_T'],v['lower_conversion_per_efold'],v['cosmic_frozen_proxy']) for k,v in q.get('best',{}).items()},flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--tradeoff',required=True);ap.add_argument('--output',required=True);ap.add_argument('--group',choices=['all','even','odd'],default='all');a=ap.parse_args()
    row=json.loads(Path(a.input).read_text())[0];q=json.loads(Path(a.tradeoff).read_text())
    scan(row,np.array(q['aa']),np.array(q['Q96']),a.output,a.group)
