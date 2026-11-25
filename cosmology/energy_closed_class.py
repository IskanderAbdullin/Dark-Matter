"""Export an energy-closed FRW spectrum without abundance rescaling; CLASS controls."""
import argparse,json,tempfile
from pathlib import Path
import numpy as np
from kinetic_psd_class import export_components,transformed_quadrature_check,transfer


def canonical(row):
    assert row['g_bath']==106.75
    s=row['final'];th=s['theta'];raw=row['frozen_cell_spectrum']
    assert raw['theta_freezeout']==th
    return dict(n=row['n'],x0=row['x0'],kappa=row['kappa'],
                final=dict(number=s['number']/th**3,omega_phi_h2_m86_g106p75=s['omega_phi_h2_m86'],
                           delta_Neff_chi2_if_stable_free_streaming=s['delta_Neff_chi2_if_stable']),
                frozen_cell_spectrum=dict(edges=(np.array(raw['edges'])/th).tolist(),f_phi=raw['f_phi']))


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',nargs='+',required=True);p.add_argument('--output',required=True)
    p.add_argument('--analytic-only',action='store_true');p.add_argument('--diagnostic',action='store_true');a=p.parse_args()
    ks=np.geomspace(.003,20.,125);zs=[0.,5.];results=[]
    with tempfile.TemporaryDirectory() as folder:
        for i,path in enumerate(a.input):
            raw=json.loads(Path(path).read_text())[0];row=canonical(raw)
            exp=export_components(row,folder,tag=str(i));alt=export_components(row,folder,tag=str(i)+'_alt',split_multiple=8.,points=8000)
            r=dict(n=row['n'],x0=row['x0'],kappa=row['kappa'],theta=raw['final']['theta'],expected_omega=row['final']['omega_phi_h2_m86_g106p75'],
                   quadrature120=transformed_quadrature_check(exp,120),quadrature240=transformed_quadrature_check(alt,240))
            if not a.analytic_only:
                r['CLASS']=transfer(row,exp,ks,zs);r['alternate_CLASS']=transfer(row,alt,ks,zs,240)
                for z in zs:
                    key=str(z);a1=r['CLASS']['ratios'][key];a2=r['alternate_CLASS']['ratios'][key];h=a1['half_power_k_per_Mpc']
                    before=ks<=(h or 10.);drift=float(np.max(np.abs(np.asarray(a2['power_ratio'])[before]/np.asarray(a1['power_ratio'])[before]-1)))
                    r.setdefault('comparisons',{})[key]=dict(half1=h,half2=a2['half_power_k_per_Mpc'],partition_drift_before_half=drift,ratio_k1=a1['P_ratio_k1'],ratio_k10=a1['P_ratio_k10'])
                    assert drift<.05
                if a.diagnostic and i==len(a.input)-1:
                    # Restrict expensive hierarchy to the suspicious enhancement region.
                    kd=np.array([.01,.05,.1,.2,.5,1.,1.5]);d={}
                    for tag,pr in [('default',{}),('late_fluid',{'ncdm_fluid_trigger_tau_over_tau_k':120.}),
                                   ('no_fluid',{'ncdm_fluid_approximation':3,'l_max_ncdm':30})]:
                        v=transfer(row,alt,kd,[0.],120,precision=pr,kmax_h=2.5)
                        d[tag]=v['ratios']['0.0']['power_ratio']
                        print('control',tag,json.dumps(d[tag]),flush=True)
                    r['fluid_diagnostic']=dict(k_per_Mpc=kd.tolist(),ratios=d)
            for e in (exp,alt):
                for c in e['components']:c.pop('file')
            r.update(export=exp,alternate_export=alt);results.append(r)
            print(json.dumps({k:v for k,v in r.items() if k not in ('CLASS','alternate_CLASS')}),flush=True)
    if len(results)==2 and not a.analytic_only:
        for z in zs:
            key=str(z);u,v=[r['CLASS']['ratios'][key] for r in results]
            h1,h2=u['half_power_k_per_Mpc'],v['half_power_k_per_Mpc']
            assert (h1 is None)==(h2 is None),(h1,h2)
            if h1 is None:
                # Both spectra must actually stay above half-power, not merely
                # lack a downward crossing (e.g. start below the threshold).
                assert min(u['power_ratio'])>.5 and min(v['power_ratio'])>.5
                pre=np.ones(len(ks),dtype=bool);shift=None
            else:
                pre=ks<=min(h1,h2);shift=abs(h1/h2-1)
            drift=float(np.max(np.abs(np.asarray(u['power_ratio'])[pre]/np.asarray(v['power_ratio'])[pre]-1)))
            print('mesh',key,h1,h2,drift,'right_censored_at',float(ks[-1]) if h1 is None else None,flush=True)
            assert (shift is None or shift<.05) and drift<.05
    Path(a.output).write_text(json.dumps(dict(k_per_Mpc=ks.tolist(),results=results),indent=2)+'\n')

if __name__=='__main__':main()
