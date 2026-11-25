"""Conditional contact-self-interaction consistency screen.

Real scalar L_int=-lambda phi^4/4!, lambda>0. Tree-level identical final
state cross section sigma(s)=lambda^2/(32*pi*s). The Hartree gap uses the
vacuum-subtracted particle variance for the ACTUAL noninteracting PSD at each
saved time. This is not a self-consistent interacting evolution, and rates
without final-state Bose factors are collision counts, not net relaxation.

Useful references:
  Arcadi et al., JHEP 08 (2019) 050, doi:10.1007/JHEP08(2019)050.
  Sakurai & Yin, JHEP 03 (2025) 202, doi:10.1007/JHEP03(2025)202.
  Du et al., PLB 845 (2023) 138161, doi:10.1016/j.physletb.2023.138161.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import brentq

GEV2_TO_CM2=0.3893793721e-27
GEV_TO_G=1.782661921e-24
MPL=1.220890e19
GYR=1e9*365.25*86400
MASS=86e-9

def sigma_over_m(lam,m=MASS):
    return lam*lam/(128*np.pi*m**3)*GEV2_TO_CM2/GEV_TO_G

def lambda_from_sigma(s,m=MASS):
    return np.sqrt(s*128*np.pi*m**3*GEV_TO_G/GEV2_TO_CM2)

def shell_nodes(edges,fp,order=8):
    l,h=np.asarray(edges[:-1]),np.asarray(edges[1:])
    u,w=roots_legendre(order)
    q=(h[:,None]+l[:,None])/2+(h-l)[:,None]*u/2
    weights=np.maximum(fp,0)[:,None]*(h-l)[:,None]*w*q*q/(4*np.pi**2)
    return q.ravel(),weights.ravel()

def variance(q,w,mu):
    return float(w@ (1/np.hypot(q,mu)))

def gap(q,w,mu0,lam):
    if lam==0:return mu0
    delta=lam/2*variance(q,w,mu0)
    if delta==0:return mu0
    # Solve for the fractional squared-mass correction. Avoid subtracting
    # two almost equal large vacuum masses at late times.
    u0=mu0*mu0
    def residual(v):
        return v-variance(q,w,np.sqrt(u0+v*delta))/variance(q,w,mu0)
    v=brentq(residual,0.,1.+1e-12,xtol=1e-14,rtol=1e-12)
    return np.sqrt(u0+v*delta)

def dilute_rate_over_H(q,w,mu,lam,M,x,total,g_bath,angular_order=12):
    # Isotropic identical relativistic elastic scattering, fixed quasiparticle
    # mass. sigma*v_M=lambda^2/(64*pi*E1*E2)*sqrt[(y-m^2)/(y+m^2)]
    # y=E1*E2-p1*p2*cos(theta). Stable at small relative momentum using rapidity.
    E=np.hypot(q,mu);a=w/E;N=w.sum()
    Href=1.66*np.sqrt(g_bath)*M*M/MPL
    Hf=np.sqrt((np.pi**2/30*g_bath)/total)
    if mu==0:
        pair=float(a.sum()**2)
    else:
        # Remove exactly empty integration nodes to limit array size.
        sel=w>0;q=q[sel];a=a[sel]
        rapid=np.arcsinh(q/mu)
        d0=2*mu*mu*np.sinh((rapid[:,None]-rapid[None,:])/2)**2
        pp=q[:,None]*q[None,:]
        u,v=roots_legendre(angular_order);pair=0.
        for c,b in zip(u,v):
            d=d0+pp*(1-c)
            pair+=float(b*.5*(a@np.sqrt(d/(d+2*mu*mu))@a))
    return float(lam*lam/(64*np.pi)*pair/max(N,1e-300)*(M*x/Href)*Hf)

def screen(row,sigmas=(.001,.01,.1,1.25)):
    specs=row.get('scalar_history_psd')
    if not specs:
        specs=[dict(x=row['final']['x'],theta=row['final']['theta'],f_phi=row['frozen_cell_spectrum']['f_phi'])]
        hist=[row['final']]
    else:
        hist=row['history'];assert len(hist)==len(specs)
    e=np.asarray(row['frozen_cell_spectrum']['edges']);out=[]
    for h,s in zip(hist,specs):
        x=s['x'];theta=s['theta'];peak=h['peak_bin_q'];mu0=row['scalar_ratio']*x
        assert abs(h['x']/x-1)<1e-12
        f=np.asarray(s['f_phi']);assert np.min(f)>-1e-7
        q,w=shell_nodes(e,f,8)
        N=w.sum();assert abs(N/h['number']-1)<1e-7
        if N<1e-10:continue
        q2,w2=shell_nodes(e,f,2);q4,w4=shell_nodes(e,f,4)
        V=variance(q,w,mu0)
        assert abs(variance(q4,w4,mu0)/V-1)<1e-4
        lam_peak=2*max(peak*peak-mu0*mu0,0)/variance(q,w,peak)
        lam_tenth=2*max((peak/10)**2-mu0*mu0,0)/variance(q,w,peak/10)
        # Exact moment inequality: <1/E> >= 1/sqrt(<q^2>+mu^2).
        jensen=N/np.hypot(h['rms_q'],mu0)
        assert V>=jensen*(1-1e-7)
        Href=1.66*np.sqrt(row['g_bath'])*row['parent_GeV']**2/MPL
        Hf=np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
        rate_coeff=V*V/N*(row['parent_GeV']*x/Href)*Hf/(64*np.pi)
        lam_rate=1/np.sqrt(rate_coeff)
        cases=[]
        for sig in sigmas:
            lam=float(lambda_from_sigma(sig))
            u=gap(q,w,mu0,lam)
            # For the conditional massive quasiparticle calculation, compare
            # independent within-cell and angular quadratures.
            r2=dilute_rate_over_H(q2,w2,u,lam,row['parent_GeV'],x,h['total_rho_Tref4'],row['g_bath'],8)
            r4=dilute_rate_over_H(q4,w4,u,lam,row['parent_GeV'],x,h['total_rho_Tref4'],row['g_bath'],16)
            cases.append(dict(sigma_over_m_cm2_g=sig,lambda_quartic=lam,mu_eff_over_Tref=u,
                gap_over_peak=float(u/peak),gap_over_bath_T=float(u/theta),
                classical_massless_rate_over_H=float(lam*lam*rate_coeff),
                conditional_massive_classical_rate_over_H=r4,
                quadrature_relative_difference=float(abs(r4/r2-1)) if r2>0 else None))
        out.append(dict(x=x,T_GeV=row['parent_GeV']*theta/x,number_over_Tref3=float(N),
            peak_q=peak,mean_q=h['mean_q'],rms_q=h['rms_q'],
            variance_over_Tref2=V,variance_Jensen_lower_bound=float(jensen),
            sigma_crit_gap_equal_peak=float(sigma_over_m(lam_peak)),
            sigma_crit_gap_tenth_peak=float(sigma_over_m(lam_tenth)),
            sigma_crit_massless_dilute_rate_H=float(sigma_over_m(lam_rate)),
            cases=cases))
    halo=[]
    for rho in [1e-25,1e-24,1e-23]:
        for v in [30.,100.]:
            for t in [1.,10.]:
                sigma=1/(rho*v*1e5*t*GYR)
                halo.append(dict(rho_g_cm3=rho,v_km_s=v,time_Gyr=t,
                    classical_one_collision_sigma_over_m_cm2_g=sigma,
                    lambda_quartic=float(lambda_from_sigma(sigma))))
    return dict(n=row['n'],kernel=row.get('kernel','original'),kappa=row['kappa'],xend=row['xend'],
        omega=row['final']['omega_phi_h2_m86'],history_available=bool(row.get('scalar_history_psd')),
        states=out,halo_classical_screen=halo,
        scope='Frozen noninteracting occupation at each sampled time; Hartree quasiparticle gap and classical elastic event rate only. No Bose factors or inverse cancellation, no mean-field stress-energy correction, no altered production, no halo relaxation proof. Contact repulsive quartic only.')

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True)
    p.add_argument('--output',required=True);a=p.parse_args()
    out=[]
    for path in a.input:
        r=json.loads(Path(path).read_text())[0]
        s=screen(r);s['input']=path;out.append(s)
        print(path,'history',s['history_available'],'min gap=peak sigma',min(x['sigma_crit_gap_equal_peak'] for x in s['states']),
            'min classical rate=H sigma',min(x['sigma_crit_massless_dilute_rate_H'] for x in s['states']),flush=True)
        Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
