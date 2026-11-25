"""Evolve actual time-integrated soft-mixing LL elastic exposure."""
import argparse,json
from pathlib import Path
from late_activation_evolve import run

def v0(row):
    h=row['final'];return h['peak_bin_q']*(3.909/row['g_bath'])**(1/3)/h['theta']
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    rows=[json.loads(Path('outputs/'+n).read_text())[0] for n in ['parent_angle240_k00685_history.json','parent_angle320_sparse_k00685.json']]
    cases=json.loads(Path('outputs/soft_mixing_loop_scan.json').read_text())
    out=[]
    for sig,z,mu in [(.1,5.,100.),(.3,8.,30.),(1.25,8.,60.)]:
        c=next(c for c in cases if c['params']['sigma']==sig and c['params']['z_on']==z and c['params']['mu_over_H_on']==mu and c['result']['loop_scheme']=='vev-curvature')
        for row in rows:
            tau=c['result']['late_z100_frozen_NR_tau']*v0(row)/v0(rows[0])
            r=run(row,z,sig,160,1e-6,30.,tau_override=tau)
            out.append(dict(params=c['params'],loop_scheme=c['result']['loop_scheme'],input_n=row['n'],result=r,
                       frozen_relative_proxy=c['result']['late_z100_frozen_N_gt2_relative_change_proxy']))
            Path(a.output).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
            print(sig,z,mu,row['n'],r['success'],r['final_moment_changes'],r['N_error'],r['kinetic_E_error'],flush=True)
if __name__=='__main__':main()
