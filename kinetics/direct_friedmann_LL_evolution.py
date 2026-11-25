"""Evolve LL contact-limit spectrum using the self-consistent Friedmann exposure."""
import argparse,json
from pathlib import Path
import numpy as np
from late_activation_evolve import run
from quartic_consistency import shell_nodes
from direct_subthreshold_activation import SCALE,VREF
from resonant_mediator_inverse import C
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mesh',type=int,required=True);p.add_argument('--input',required=True);p.add_argument('--exposure',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    row=json.loads(Path(a.input).read_text())[0];v=next(v for v in json.loads(Path(a.exposure).read_text()) if v['mesh']==a.mesh)
    r=run(row,v['params']['z_on'],v['sigma0'],160,1e-6,100.,tau_override=v['LL_tau'])
    q,w=shell_nodes(np.array(row['frozen_cell_spectrum']['edges']),np.maximum(np.array(row['frozen_cell_spectrum']['f_phi']),0.),8)
    q=q/row['final']['peak_bin_q'];v0=v['full_friedmann']['v0']
    # 1% amplitude-squared correction at r=101 and r=10: if both
    # momenta are within qmax, vrel<=2 qmax*v0*r, so compare to detuning.
    tails=[]
    for R in [101.,10.,5.]:
        qmax=np.sqrt(v['a']*(1/np.sqrt(.99)-1))*VREF/(2*v0*R)
        tails.append(dict(r=R,qmax_over_peak=float(qmax),number_tail=float(sum(w[q>qmax])/sum(w)),
                        kinetic_tail=float(sum(w[q>qmax]*q[q>qmax]**2)/sum(w*q*q))))
    Path(a.output).write_text(json.dumps(dict(mesh=a.mesh,exposure=v,kinetic=r,contact_validity=tails),indent=2,allow_nan=False)+'\n')
    print(a.mesh,r['success'],r['final_moment_changes'],r['N_error'],r['kinetic_E_error'],tails,flush=True)
