"""Scan explicit one-loop eta potential matching choices in soft mixing."""
import argparse,json
from pathlib import Path
from soft_mixing_activation import scales,evolve

def scan(row,output,refined_bath=True):
    out=[]
    for sigma,z,muH in [(.1,5.,30.),(.1,5.,100.),(.1,8.,30.),(.1,8.,100.),(.3,8.,30.),(.3,8.,100.),(1.25,8.,20.),(1.25,8.,60.),(1.25,10.,30.),(1.25,10.,60.)]:
        p0=scales(sigma,z_on=z,dneff=.1,muH=muH)
        if p0['overlap_ratio']<1.05:continue
        p=scales(sigma,z_on=z,dneff=.1,muH=muH,X=1.05*max(p0['X_loop_floor'],p0['X_stability_floor']))
        for scheme in ['origin-vev','vev-curvature']:
            r=evolve(row,p,1e-8,loop_scheme=scheme,refined_bath=refined_bath)
            out.append(dict(params=p,result=r))
            Path(output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print(sigma,z,muH,scheme,r['first_sigma_fraction_crossings']['0.5'],r['today_sigma_LL'],r['conditional_local_LL_rayleigh_efolds'],r['late_z100_frozen_N_gt2_relative_change_proxy'],r['eta_CW_over_X'],flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--legacy-bath',action='store_true');a=p.parse_args()
    scan(json.loads(Path(a.input).read_text())[0],a.output,not a.legacy_bath)
