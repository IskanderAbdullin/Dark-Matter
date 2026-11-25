"""Targeted hierarchy convergence around the physical half-power crossing."""
import argparse,json,tempfile
from pathlib import Path
import numpy as np
from energy_closed_class import canonical
from kinetic_psd_class import export_components,transfer


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    row=canonical(json.loads(Path(a.input).read_text())[0])
    ks=np.unique(np.r_[.003,.01,.03,np.geomspace(.1,1.5,65)])
    with tempfile.TemporaryDirectory() as folder:
        e=export_components(row,folder,split_multiple=8.,points=8000)
        cases=[('default',120,{}),('full_l30',120,{'ncdm_fluid_approximation':3,'l_max_ncdm':30}),
               ('full_l45',240,{'ncdm_fluid_approximation':3,'l_max_ncdm':45,'tol_perturbations_integration':3e-6})]
        out={}
        for label,bins,precision in cases:
            v=transfer(row,e,ks,[0.,5.],bins,precision=precision,kmax_h=2.5)
            for z,s in v['ratios'].items():
                # Transfer helper's k=10 interpolation is outside this run's range.
                s.pop('P_ratio_k10',None)
                print(label,z,'half',s['half_power_k_per_Mpc'],'P1',s['P_ratio_k1'],flush=True)
            out[label]=v
        checks={}
        for z in ['0.0','5.0']:
            a1=out['full_l30']['ratios'][z];a2=out['full_l45']['ratios'][z]
            h1,h2=a1['half_power_k_per_Mpc'],a2['half_power_k_per_Mpc'];assert h1 and h2
            r1=np.asarray(a1['power_ratio']);r2=np.asarray(a2['power_ratio'])
            before=ks<=min(h1,h2)
            checks[z]=dict(half_shift=abs(h1/h2-1),max_pre_half_drift=float(np.max(np.abs(r1[before]/r2[before]-1))),
                           max_all_drift=float(np.max(np.abs(r1/r2-1))))
        Path(a.output).write_text(json.dumps(dict(k_per_Mpc=ks.tolist(),runs=out,checks=checks),indent=2)+'\n')
        print(json.dumps(checks),flush=True)
        assert all(c['half_shift']<.03 and c['max_pre_half_drift']<.03 for c in checks.values()),checks
