"""Direct mediator: self-consistent homogeneous Friedmann test.

Use exact decoupled eta energy/pressure tabulation; solve chi with the
algebraic Friedmann constraint (kinetic term depends on H). Standard
matter/radiation/Lambda kept fixed; no retuning H0 or cosmological fit.
The very small late A/S medium contribution is omitted here.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import solve_ivp,trapezoid,cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from direct_subthreshold_activation import VREF,SCALE,RHO,sunset_threshold
from direct_subthreshold_background import energy
from direct_activation_consistency import rho_tot,MPL_REDUCED
from soft_mixing_activation import variance_table
from late_activation_overlap import H,H0,Om,Or
from resonant_mediator_inverse import HBAR,M,C,SEC_GYR
from quartic_consistency import lambda_from_sigma
from velocity_only_exposure import RHO0

def eta_table(order=48):
    u,w=roots_legendre(order);e=np.array([0.,1e-5,1e-4,.001,.01,.1,1.,5.,15.,50.])
    q=((e[:-1,None]+e[1:,None])/2+np.diff(e)[:,None]*u/2).ravel()
    f=q*q*(np.diff(e)[:,None]*w/2).ravel()/np.expm1(q)/(2*np.pi**2)
    j=np.geomspace(1e-8,1e5,1200);E=np.hypot(j[:,None],q[None,:])
    vals=[E@f,(1/E)@(f*q*q)/3,(1/E)@f]
    ips=[PchipInterpolator(np.log(j),np.log(v)) for v in vals]
    def get(y):
        y=np.asarray(y);t=np.log(np.clip(y,j[0],j[-1]))
        return [np.exp(p(t)) for p in ips]
    return get

def solve_friedmann(row,p,rtol=1e-9,n=24000,order=48,return_history=False):
    et=eta_table(order);X=p['X'];v=p['vev_eV'];mu=p['mu_eV'];cw=p['eta_CW_over_X']
    td0=p['Td_on_eV']/(p['z_on']+1);m0=p['m_eta_at_vev_eV'];N=p['Neta']
    rho0=3*MPL_REDUCED**2*(H0*HBAR)**2
    def local(t,x,dx):
        r=np.exp(-t);td=td0*r;ar,ap,iv=et(m0*abs(x)/td)
        er=N*td**4*ar;ep=N*td**4*ap
        xx=x*x;lg=np.log(max(xx,1e-300))
        U=X*(.25*(xx-1)**2+cw*(xx*xx*lg-1.5*xx*xx+2*xx-.5))
        dU=X*x*(xx-1+4*cw*(1+xx*(lg-1)))+N*td*td*m0*m0*x*iv
        rb=rho0*(Om*r**3+Or*r**4+1-Om-Or)
        den=3*MPL_REDUCED**2-.5*v*v*dx*dx
        assert den>0.
        hh=np.sqrt((rb+U+er)/den)
        K=.5*v*v*hh*hh*dx*dx
        dh=-.5*(rho0*(Om*r**3+4/3*Or*r**4)+2*K+er+ep)/(MPL_REDUCED**2*hh*hh)
        return hh,dh,dU,K+U+er,K-U+ep
    def fun(t,y):
        x,dx=y;hh,dh,dU,_,_=local(t,x,dx)
        return [dx,-(3+dh)*dx-dU/(v*v*hh*hh)]
    t0=-np.log(1e14);ts=np.r_[np.linspace(t0,-np.log(101.),500),np.linspace(-np.log(101.),0.,n)[1:]]
    sol=solve_ivp(fun,(t0,0.),[1e-8,0.],t_eval=ts,rtol=rtol,atol=rtol*1e-10,method='DOP853',max_step=.04,first_step=.01)
    assert sol.success,sol.message
    r=np.exp(-ts);x,dx=sol.y
    loc=np.array([local(t,xx,dd) for t,xx,dd in zip(ts,x,dx)]);hh=loc[:,0]/HBAR;rho=loc[:,3];pr=loc[:,4]
    l=r<=101;h=r<=3;A4=trapezoid(x[h]**4/hh[h],ts[h])
    cosmic=trapezoid(r[l]**4*x[l]**4/hh[l],ts[l])
    _,Tr=variance_table(row);v0=row['final']['peak_bin_q']*Tr/M*C
    R=rho[l]/r[l]**4
    target=R[0]+cumulative_trapezoid((rho[l]-3*pr[l])/r[l]**4,ts[l],initial=0.)
    summary=dict(active_time_Gyr=float(A4/SEC_GYR),A4=float(A4),cosmic=float(cosmic),v0=float(v0),
      max_extra_fraction_zle100=float(max(rho[l]/(3*MPL_REDUCED**2*(hh[l]*HBAR)**2))),
      max_relative_H_change_zle100=float(max(abs(hh[l]/H(r[l]-1)-1))),H0_ratio=float(hh[-1]/H0),
      continuity_residual=float(max(abs(R-target))/max(abs(R))),max_field=float(max(abs(x))),today_field=float(x[-1]),
      z_half=float(r[np.flatnonzero(abs(x)**4>=.5)[0]]-1),rtol=rtol,n=n,order=order)
    if return_history:
        return summary, dict(log_a=ts, r=r, x=x, dx=dx, H=hh, rho_extra=rho, pressure_extra=pr)
    return summary

def evaluate(row,s,b,mesh):
    p=s['params'];rr=solve_friedmann(row,p,n=12000 if mesh==240 else 20000,rtol=1e-8 if mesh==240 else 2e-9,order=32 if mesh==240 else 48)
    sig=1/(RHO*VREF*b['Q']*rr['A4']);g2=float(lambda_from_sigma(sig))*b['a']*SCALE
    tau=sig*RHO0*rr['v0']*rr['cosmic']
    return dict(mesh=mesh,params=p,a=b['a'],Q=b['Q'],sigma0=float(sig),g2_eV2=float(g2),
        LL_tau=float(tau),cosmic_proxy=float(5.756e7*tau),
        max_loop_over_X_scales=max(2*abs(sunset_threshold(g2,q*M)['total_eV4'])/p['X'] for q in [.5,1,2]),
        full_friedmann=rr,scope='Homogeneous self-consistent expansion, renormalized trigger, no perturbations, no H0 refit or primordial S.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--input2',required=True);ap.add_argument('--scans',nargs='+',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    rows=[json.loads(Path(f).read_text())[0] for f in [a.input,a.input2]]
    scans=sum([json.loads(Path(f).read_text()) for f in a.scans],[])
    # Candidate at lowest radiative ratio under exact energy <=2%.
    candidates=[(s,b) for s in scans for b in s['screen'].get('rows',[]) if b['exact_chi_eta_fraction']<.02 and b['cosmic_proxy']<.01 and b['LL_sigma_m_1000kms']<.2]
    s,b=min(candidates,key=lambda z:z[1]['max_loop_over_X_scales'])
    out=[]
    for row,mesh in zip(rows,[240,320]):
        v=evaluate(row,s,b,mesh);out.append(v)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
        print('FRIEDMANN',mesh,v['params']['z_on'],v['params']['dneff'],v['a'],v['cosmic_proxy'],v['max_loop_over_X_scales'],v['full_friedmann'],flush=True)
