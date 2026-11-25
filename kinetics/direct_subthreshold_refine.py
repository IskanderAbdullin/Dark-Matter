"""Refined single-light-state direct mediator test and early diagnostics."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import qmc
from direct_subthreshold_activation import solve,diagnose,sunset_threshold,direct_recombination,VREF,SCALE,RHO
from direct_activation_consistency import CGS_RHO_TO_EV4
from subthreshold_recombination import radial_bose,avg_rate
from soft_mixing_activation import variance_table
from early_bose_bootstrap import rearrange,N0
from quartic_consistency import shell_nodes
from resonant_mediator_inverse import M,C,HBAR,SEC_GYR
from late_activation_overlap import H
from late_activation_evolve import run as kinetic_run

def early(row,rr,b,p,order=2):
    """Conditional early collision bounds and frozen forward self-energy.

    Start at post-portal r=1e14; A number <=2*Nref, A cap=2.1 times
    frozen maximum. S,chi occupations assumed <=1 for direct number
    production rates (not proved for the two species here).
    """
    r,t,x=rr['r'],rr['t'],rr['x'];early=r>=1001
    r,t,x=r[early],t[early],x[early];h=.5*HBAR*H(r-1)
    n=N0*r**3;F=2.1*max(row['frozen_cell_spectrum']['f_phi'])
    det=b['a']*SCALE;ms2=4*M*M-det;ms=np.sqrt(ms2)
    g2=b['g2_eV2']*x*x;v=p['vev_eV']
    # s >=4M^2, t,u<=0 in elastic scattering. Triangle bound.
    amp=g2*(1/det+5/ms2)
    I,V=rearrange(2*n,F,M)
    mode=amp**2*(1+2*F)*I*I/(64*np.pi**3*M*h)
    l1=4*amp**2*(1+2*F)*V*V/(32*np.pi*n*h)
    # For AA -> SS, t,u<=0; initial S cap of 1 assumed.
    ass=2*g2/(M*M)
    ss=2*ass**2*4*V*V/(32*np.pi*n*h)
    # AA -> chi S via direct dimensionless kappa. Outgoing caps <=1.
    kapp2=b['g2_eV2']/v**2
    chi=4*kapp2*V*V/(32*np.pi*n*h)
    # 3A -> A S: all tree denominators bounded away from poles at
    # fixed masses because mS < 2M. For the mixed-sign external
    # pair S-A the closest virtual A denominator is M^2-(mS-M)^2.
    # The two A incoming pair denominator >= det, mixed A pair <=0.
    ad=M*M-(ms-M)**2
    bound5=4*np.sqrt(g2)*(3*g2/ms2+3*g2/det)/ad
    # Two-body integrated phase space <=1/(8pi), initial 3! and 2^3.
    rec=bound5**2*V**3*(1+F)*2/(384*np.pi*n*h)
    def integrate(y):return float(trapezoid(y,t))
    # Actual frozen PSD forward-scattering estimate. Use common-energy
    # free dispersion and leading-density real self-energies; no
    # self-consistent quasiparticle solution. Sample all occupied q and p=0.
    q,w=shell_nodes(np.asarray(row['frozen_cell_spectrum']['edges']),np.asarray(row['frozen_cell_spectrum']['f_phi']),order)
    good=w>0;q,w=q[good],w[good];w/=sum(w)
    _,Tr=variance_table(row)
    data=[]
    for R,X in zip(r[::5],x[::5]):
        P=q*Tr*R;E=np.hypot(M,P);Y=np.arcsinh(P/M);g=b['g2_eV2']*X*X;Nr=N0*R**3
        ext=np.r_[0.,P];ye=np.arcsinh(ext/M)
        A=det+4*M*M*np.sinh((ye[:,None]-Y[None,:])/2)**2
        B=4*ext[:,None]*P[None,:]
        avg=np.log1p(B/A)/np.where(B>0,B,1.)
        avg[0]=1/A[0]
        # Forward t=0; u<=0; use absolute triangle estimate with
        # generous unit multiplicity, no cancellations assumed.
        VA=Nr*np.sum(w/E)
        PiA=g*(Nr*np.max(avg@(w/E))+5*VA/ms2)
        PiS=2*g*Nr*np.sum(w/(E*(det+4*P*P)))
        data.append(dict(r=float(R),x=float(X),PiA_bound_over_det=float(PiA/det),
                         PiS_rest_bound_over_det=float(PiS/det),VA_eV2=float(VA)))
    res=dict(start_r=float(r[0]),end_r=float(r[-1]),FA_cap=float(F),initial_F_max=float(F/2.1),
       max_early_abs_x=float(max(abs(x))),max_elastic_mode_gain_over_H=float(max(mode)),
       integrated_elastic_mode_gain_bound=integrate(mode),elastic_L1_over_Nref_bound=integrate(l1),
       integrated_SS_number_upper_if_FS_le1=integrate(ss),
       integrated_chiS_event_upper_if_caps_le1=integrate(chi),
       integrated_3A_to_AS_event_upper_if_FS_le1=integrate(rec),
       max_forward_PiA_over_det=max(x['PiA_bound_over_det'] for x in data),
       max_forward_PiS_rest_over_det=max(x['PiS_rest_bound_over_det'] for x in data),
       mass_diagnostics=data,
       scope='Postportal conditional A occupancy/number caps, stable vacuum poles, S/chi final caps for production only. No S mode bootstrap, primordial S/chi spectrum, or proof of the nonlocal medium approximation.')
    return res

def halo_rec(b,powers=(16,18)):
    a=b['a'];g=np.sqrt(b['g2_eV2']);ms=np.sqrt(4*M*M-a*SCALE);lam=3*g*g/ms**2
    alpha=1e-6;scale=M*VREF/C;n=RHO*CGS_RHO_TO_EV4/M
    xs,cdf,meta=radial_bose(alpha)
    radial=lambda v:np.interp(v,cdf,xs)
    f=lambda q:1/np.expm1(alpha+(q/scale)**2)
    out=[]
    for power in powers:
        for seed in [137,138]:
            s=qmc.Sobol(d=8,scramble=True,seed=seed).random_base2(power)
            R=avg_rate(s,radial,scale,1.,g,lam,M,M+1.,ms,f,False)
            out.append(dict(power=power,seed=seed,initial_S_event_fraction_10Gyr=R['coefficient']*n*n*10*SEC_GYR/HBAR,rate=R))
    return dict(radial_checks=meta,runs=out)

def pick(scans):
    # Lower energy-fraction candidate preferred over a marginal 4% case.
    # Pick best MS-M threshold at energy fraction <1% and <4%.
    all=[]
    for s in scans:
        for b in s['screen'].get('rows',[]):
            if b['cosmic_proxy']<.0095 and b['LL_sigma_m_1000kms']<.2:
                all.append((s,b))
    for lim in [.01,.04]:
        pool=[(s,b) for s,b in all if b['max_tree_plus_eta_fraction_zle100']<lim]
        if pool:
            yield min(pool,key=lambda z:z[1]['max_loop_over_X_scales'])

def refine(row,s,b,mesh,do_kinetic=True):
    p=s['params'];scheme=s['scheme']
    rr=solve(row,p,scheme,g2=b['g2_eV2'],n=22000,rtol=2e-9)
    qq=np.array([b['Q']]);aa=np.array([b['a']])
    d=diagnose(row,p,rr,aa,qq);v=d['rows'][0]
    # One extra fixed-point update because g2 depends on roll integral.
    rr=solve(row,p,scheme,g2=v['g2_eV2'],n=22000,rtol=2e-9)
    d=diagnose(row,p,rr,aa,qq);v=d['rows'][0]
    det=v['a']*SCALE;g2=v['g2_eV2'];ms2=4*M*M-det
    cutoff=(4*p['X']*ms2*(16*np.pi**2)**2/(3*g2))**.25
    out=dict(params=p,scheme=scheme,mesh=mesh,screen=v,z_half=d['z_half'],active_time_Gyr=d['active_time_Gyr'],
       g_eV=float(np.sqrt(g2)),mS_eV=float(np.sqrt(ms2)),binding_eV=float(det/(2*M+np.sqrt(ms2))),
       relative_mass_detuning=float(det/(4*M*M)),wilsonian_figure8_cutoff_eV=float(cutoff),
       one_loop_mass_over_det=float(g2/(16*np.pi**2*det)),kappa=float(np.sqrt(g2)/p['vev_eV']),
       ms_thresholds=[sunset_threshold(g2,Q*M) for Q in [.5,1.,2.]],
       roll={k:v for k,v in rr.items() if k not in ['t','r','x','dx']},
       early=early(row,rr,v,p),cosmic_recombination=direct_recombination(row,p,rr,v))
    if mesh==240:out['halo_recombination']=halo_rec(v)
    if do_kinetic:
        out['kinetic']=kinetic_run(row,p['z_on'],v['sigma0'],160,1e-6,90.,tau_override=v['LL_tau'])
    print('DIRECT REFINE',mesh,p['z_on'],p['dneff'],p['Xfrac'],p['mu_over_H_on'],scheme,
          {k:v[k] for k in ['a','sigma0','cosmic_proxy','loop_over_X_at_M','max_loop_over_X_scales','max_tree_plus_eta_fraction_zle100']},
          {k:out['early'][k] for k in ['integrated_elastic_mode_gain_bound','elastic_L1_over_Nref_bound','integrated_3A_to_AS_event_upper_if_FS_le1','max_forward_PiA_over_det','max_forward_PiS_rest_over_det']},
          out.get('kinetic',{}).get('final_moment_changes'),flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--input2',required=True);ap.add_argument('--scans',nargs='+',required=True);ap.add_argument('--output',required=True);ap.add_argument('--index',type=int,default=-1);a=ap.parse_args()
    rows=[json.loads(Path(f).read_text())[0] for f in [a.input,a.input2]]
    scans=sum([json.loads(Path(f).read_text()) for f in a.scans],[])
    out=[]
    for i,(s,b) in enumerate(pick(scans)):
        if a.index>=0 and i!=a.index:continue
        runs=[]
        for mesh,row in zip([240,320],rows):
            runs.append(refine(row,s,b,mesh))
            Path(a.output).write_text(json.dumps(out+[dict(index=i,runs=runs)],indent=2,allow_nan=False)+'\n')
        out.append(dict(index=i,runs=runs))
