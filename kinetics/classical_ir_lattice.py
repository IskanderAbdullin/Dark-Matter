"""Classical-statistical IR-only real scalar pilot from the actual PSD.

Gaussian random initial conditions in periodic box, with mean-field
frequency fixed self-consistently from the discrete ensemble covariance.
Only highly occupied modes below an explicit physical momentum cut are
represented; UV/low-occupancy gas and interactions with it are absent.
This is not a complete relativistic scalar relic evolution.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from coherent_lattice_trial import accel,advance,lattice_k2
from radial_lattice_map import occupation_radial
from quartic_consistency import shell_nodes,lambda_from_sigma,gap,variance,MPL

def run(row,sigma,n,L,dt,end,cut,seed,fmin=100.,mapping="direct"):
    e=np.array(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    h=row['final'];x=h['x'];p0=h['peak_bin_q'];m0=1e-12*x
    lam=float(lambda_from_sigma(sigma));alpha=(m0/p0)**2;beta=1-alpha
    dx=L/n;k2=lattice_k2(n,L);k=np.sqrt(k2)
    q=p0*k
    ix=np.searchsorted(e,q,side='right')-1
    fg=np.where((ix>=0)&(ix<len(f))&(k<=cut),f[np.clip(ix,0,len(f)-1)],0.)
    fg=np.where(fg>=fmin,fg,0.);fg[0,0,0]=0.
    mapinfo={}
    if mapping=="radial":
        fg,mapinfo=occupation_radial(k2,e,f,p0,L,cut,fmin)
    cf=lam/(6*beta*dx**3)
    # Expected Hartree variance in this finite lattice ensemble.
    def resid(s):
        om=np.sqrt(k2+s)
        return s-alpha-3*beta*cf*np.mean(fg/om)
    upper=alpha+3*beta*cf*np.mean(fg/np.sqrt(k2+alpha))
    m2=brentq(resid,alpha,upper*1.0001,xtol=1e-15,rtol=1e-12)
    om=np.sqrt(k2+m2)
    rng=np.random.default_rng(seed)
    def gaussian(variance):
        z=np.fft.fftn(rng.standard_normal(k.shape),norm='ortho')*np.sqrt(variance)
        return np.fft.ifftn(z,norm='ortho').real
    y=gaussian(cf*fg/om);vel=gaussian(cf*fg*om)
    a=np.empty_like(y);accel(y,a,1/dx**2,alpha,beta);v=vel-dt*a/2
    assert dt*np.sqrt(k2.max()+alpha+3*beta*np.max(y*y))<1.5
    Nv=p0**3*6*beta/lam;Uv=p0**4*6*beta/lam
    # Continuum cut and representable number/energy quantified independently.
    qc,wc=shell_nodes(e,f,16)
    ff=f[np.clip(np.searchsorted(e,qc)-1,0,len(f)-1)]
    sel=(qc<=cut*p0)&(ff>=fmin)
    qall=qc;wall=wc
    mfull=gap(qall,wall,m0,lam);Ef=np.hypot(qall,mfull)
    mcut=gap(qc[sel],wc[sel],m0,lam);Ecut=np.hypot(qc[sel],mcut)
    Ufull=float(wall@Ef-lam*variance(qall,wall,mfull)**2/8)
    Ucut=float(wc[sel]@Ecut-lam*variance(qc[sel],wc[sel],mcut)**2/8)
    Uexp=p0**4*(float(np.sum(fg*om)/L**3)-lam*(float(np.sum(fg/om)/L**3))**2/8)
    Nexp=p0**3*float(np.sum(fg)/L**3)
    bins=np.r_[0.,np.geomspace(.1,max(12.,np.sqrt(k2.max())*1.01),29)]
    bidx=np.clip(np.searchsorted(bins,k,side="right")-1,0,len(bins)-2)
    counts=np.bincount(bidx.ravel(),minlength=len(bins)-1)
    def stats(t):
        accel(y,a,1/dx**2,alpha,beta);u=v+dt*a/2
        kin=.5*np.mean(u*u)
        grad=.5/dx**2*sum(np.mean((np.roll(y,1,axis=ax)-y)**2) for ax in range(3))
        pot=np.mean(alpha*y*y/2+beta*y**4/4)
        mH2=alpha+3*beta*np.mean(y*y)
        omega=np.sqrt(k2+mH2)
        z=np.fft.fftn(y,norm='ortho');p=np.fft.fftn(u,norm='ortho')
        Nmode=(abs(p)**2+omega**2*abs(z)**2)/(2*omega*y.size)*Nv
        Eqmode=Nmode*omega*p0
        N=float(Nmode.sum());Q2=float(np.sum(Nmode*q*q)/N)
        bN=np.bincount(bidx.ravel(),weights=Nmode.ravel(),minlength=len(bins)-1)
        bU=np.bincount(bidx.ravel(),weights=Eqmode.ravel(),minlength=len(bins)-1)
        return dict(radial_number_proxy=bN.tolist(),radial_energy_proxy=bU.tolist(),t=t,U=float((kin+grad+pot)*Uv),N_hartree_proxy=N,
           q2_proxy=Q2,mH_over_peak=float(np.sqrt(mH2)),
           zero_N_proxy_fraction=float(Nmode[0,0,0]/N),
           low_k_N_proxy_fraction=float(Nmode[(k>0)&(k<.5)].sum()/N),
           beyond_initial_cut_N_proxy_fraction=float(Nmode[k>cut].sum()/N),
           UV_energy_proxy_fraction=float(Eqmode[k>.75*np.sqrt(k2.max())].sum()/Eqmode.sum()),
           mean_field=float(y.mean()),field_variance=float(np.var(y)))
    hist=[];start=time.monotonic();steps=int(round(end/dt));every=max(1,int(round(5/dt)))
    for i in range(0,steps+1,every):
        hist.append(stats(i*dt))
        if i>=steps:break
        advance(y,v,a,1/dx**2,alpha,beta,dt,min(every,steps-i))
    if hist[-1]['t']<steps*dt:hist.append(stats(steps*dt))
    U0=hist[0]['U'];N0=hist[0]['N_hartree_proxy']
    TH=MPL*x/(1.66*np.sqrt(row['g_bath'])*row['parent_GeV'])*np.sqrt(np.pi**2/30*row['g_bath']/h['total_rho_Tref4'])
    out=dict(radial_proxy_bins_over_peak=bins.tolist(),radial_proxy_counts=counts.tolist(),mapping=mapping,mapping_info=mapinfo,sigma=sigma,lambda_quartic=lam,n=n,L=L,dt=dt,cut_over_peak=cut,seed=seed,fmin=fmin,
        peak_q=p0,m_eff_expected_over_peak=np.sqrt(m2),N_continuum_full=float(wall.sum()),
        U_continuum_full=Ufull,N_continuum_cut=float(wc[sel].sum()),U_continuum_cut=Ucut,
        N_lattice_expected=Nexp,U_lattice_expected=Uexp,
        lattice_expected_N_over_cut=Nexp/wc[sel].sum(),lattice_expected_U_over_cut=Uexp/Ucut,
        actual_initial_N_over_expected=N0/Nexp,actual_initial_U_over_expected=U0/Uexp,
        lattice_modes=int((fg>0).sum()),maximum_energy_relative_drift=max(abs(h['U']/U0-1) for h in hist),
        final=hist[-1],history=hist,wall_seconds=time.monotonic()-start,
        final_H_time=end/(p0*TH),scope='IR only classical statistical pilot; quasiparticle number proxy, not exact charge; omitted UV tail and production, no claim of relic abundance or condensation')
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--sigma',type=float,default=.001);p.add_argument('--n',type=int,default=64)
    p.add_argument('--L',type=float,default=24.);p.add_argument('--dt',type=float,default=.025)
    p.add_argument('--end',type=float,default=180.);p.add_argument('--cut',type=float,default=6.)
    p.add_argument('--seed',type=int,default=42);p.add_argument('--mapping',choices=['radial','direct'],default='radial');a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0]
    r=run(row,a.sigma,a.n,a.L,a.dt,a.end,a.cut,a.seed,mapping=a.mapping)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k!='history'}),flush=True)
