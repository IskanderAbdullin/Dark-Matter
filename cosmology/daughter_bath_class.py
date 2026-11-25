"""Conditional transfer/export after the massive daughter returns heat to the SM bath.

No artificial PSD normalization. Reuse the synchronous-gauge CLASS
pipeline, and independently compute a trajectory free-streaming diagnostic.
The latter is NOT a matter transfer function or an observational likelihood.
"""
import argparse
import json
import tempfile
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre
from kinetic_psd_class import transformed_quadrature_check,transfer
from three_scale_psd import export_three


def canonical(row):
    assert row['g_bath']==106.75 and row['xend']>=60.
    raw=row['frozen_cell_spectrum'];s=row['final'];theta=s['theta']
    assert raw['theta_freezeout']==theta
    assert abs(s['rho_chi2_Tref4'])/s['total_rho_Tref4']<1e-8
    assert abs(s['mother_rho_Tref4'])/s['total_rho_Tref4']<1e-8
    previous=row['history'][-2]
    omega_drift=abs(s['omega_phi_h2_m86']/previous['omega_phi_h2_m86']-1)
    assert omega_drift<.002,omega_drift
    # Under this model all chi2 energy is deposited in the SM thermal bath
    # BEFORE neutrino decoupling; theta already accounts for entropy dilution.
    # No surviving extra radiation species from chi2 is included.
    T_end_GeV=theta*row['parent_GeV']/row['xend']
    assert T_end_GeV>0.1,T_end_GeV
    return dict(n=row['n'],x0=row['x0'],kappa=row['kappa'],
                final=dict(number=s['number']/theta**3,omega_phi_h2_m86_g106p75=s['omega_phi_h2_m86'],
                           delta_Neff_chi2_if_stable_free_streaming=0.),
                frozen_cell_spectrum=dict(edges=(np.array(raw['edges'])/theta).tolist(),f_phi=raw['f_phi']),
                freezeout=dict(T_bath_GeV=T_end_GeV,theta=theta,omega_drift=omega_drift,
                               residual_chi2_fraction=s['matter_fraction_chi2']))


def shell_nodes(row,order=8):
    e=np.asarray(row['frozen_cell_spectrum']['edges']);f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    x,w=roots_legendre(order);l,h=e[:-1,None],e[1:,None]
    q=l+(h-l)*(x+1)/2
    weights=f[:,None]*q*q*(h-l)*w/(4*np.pi**2)
    assert abs(weights.sum()/row['final']['number']-1)<1e-8
    return q.ravel(),weights.ravel()


def free_streaming(row,order=8,time_order=160):
    q,w=shell_nodes(row,order);w=w/w.sum()
    h=.674;omega_b=.02237;omega_dm=row['final']['omega_phi_h2_m86_g106p75']
    Om_b=omega_b/h**2;Om_dm=omega_dm/h**2
    Om_r=2.4728e-5/h**2*(1+.22710731766*3.044)
    Om_l=1-Om_b-Om_dm-Om_r
    T0=2.7255*8.617333262e-5;tr=(3.909/106.75)**(1/3)
    a0=T0*tr/(row['freezeout']['T_bath_GeV']*1e9)
    aeq=Om_r/(Om_dm+Om_b)  # equal-density CDM reference equality
    anr=q*T0*tr/86.
    norm=float(w@np.hypot(1.,anr))
    def distance(aend):
        x,g=roots_legendre(time_order)
        lo,hi=np.log(a0),np.log(aend);a=np.exp((lo+hi)/2+(hi-lo)*x/2)
        # a^4 E(a)^2 from exact massive-particle background energy moment.
        dm=(np.hypot(a[:,None],anr[None,:])@w)*Om_dm/norm
        denom=np.sqrt(Om_r+Om_b*a+dm+Om_l*a**4)
        v=anr[None,:]/np.hypot(a[:,None],anr[None,:])
        return (299792.458/(100*h))*((g*(hi-lo)/2*a/denom)@v)
    def stats(L):
        idx=np.argsort(L);cdf=np.cumsum(w[idx])
        return dict(mean_Mpc=float(w@L),rms_Mpc=float(np.sqrt(w@(L*L))),
                    quantiles_Mpc={str(p):float(np.interp(p,cdf,L[idx])) for p in [.5,.9,.99,.999]},
                    fractions_above={str(t):float(w[L>t].sum()) for t in [.01,.06,.1,1.,10.]})
    L_eq=distance(aeq);L_today=distance(1.)
    return dict(a_freezeout=a0,a_eq_CDM_reference=aeq,through_equality=stats(L_eq),through_today=stats(L_today),
                number_weighted_speed_today_km_s=float(299792.458*(w@(anr/np.hypot(1.,anr)))),
                scope='homogeneous collisionless trajectories from freezeout; not a linear transfer bound')


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    p.add_argument('--analytic-only',action='store_true');p.add_argument('--full-hierarchy',action='store_true')
    p.add_argument('--kmax',type=float,default=20.);a=p.parse_args()
    ks=np.geomspace(.003,a.kmax,125);zs=[0.,5.];results=[]
    with tempfile.TemporaryDirectory() as folder:
        for i,path in enumerate(a.input):
            raw=json.loads(Path(path).read_text())[0];row=canonical(raw)
            fs=free_streaming(row);fs2=free_streaming(row,12,240)
            drift=max(abs(fs[key][stat]/fs2[key][stat]-1) for key in ['through_equality','through_today'] for stat in ['mean_Mpc','rms_Mpc'])
            assert drift<1e-5,drift
            result=dict(input=path,n=row['n'],kappa=row['kappa'],expected_omega=row['final']['omega_phi_h2_m86_g106p75'],
                        freezeout=row['freezeout'],free_streaming=fs,free_streaming_moment_drift=drift)
            print('free_streaming',json.dumps(result),flush=True)
            # Positive three-tracer export, independently varied partition.
            e=export_three(row,folder,tag=str(i))
            e2=export_three(row,folder,tag=str(i)+'_alt',cold_multiple=10.,middle_factor=.8,points=9000)
            result.update(export=e,alternate_export=e2,quadrature120=transformed_quadrature_check(e,120),
                          quadrature240=transformed_quadrature_check(e2,240))
            if not a.analytic_only:
                precision={'ncdm_fluid_approximation':3,'l_max_ncdm':30} if a.full_hierarchy else None
                precision2={'ncdm_fluid_approximation':3,'l_max_ncdm':45,'tol_perturbations_integration':3e-6} if a.full_hierarchy else None
                result['CLASS']=transfer(row,e,ks,zs,precision=precision,kmax_h=max(2.5,a.kmax/.674*1.15))
                result['alternate_CLASS']=transfer(row,e2,ks,zs,240,precision=precision2,kmax_h=max(2.5,a.kmax/.674*1.15))
                comparisons={}
                for z in zs:
                    a1=result['CLASS']['ratios'][str(z)];a2=result['alternate_CLASS']['ratios'][str(z)]
                    v1=np.array(a1['power_ratio']);v2=np.array(a2['power_ratio'])
                    comparisons[str(z)]=dict(half1=a1['half_power_k_per_Mpc'],half2=a2['half_power_k_per_Mpc'],
                        minimum_ratio=float(min(v1)),full_range_drift=float(np.max(np.abs(v1/v2-1))),P1=a1['P_ratio_k1'],
                        P10=a1['P_ratio_k10'] if a.kmax>=10 else None)
                result['comparisons']=comparisons
            for ex in (e,e2):
                for c in ex['components']:c.pop('file')
            results.append(result)
            Path(a.output).write_text(json.dumps(dict(k_per_Mpc=ks.tolist(),results=results),indent=2)+'\n')
            print('export',json.dumps({k:v for k,v in result.items() if k not in ('CLASS','alternate_CLASS','free_streaming')}),flush=True)


if __name__=='__main__':main()
