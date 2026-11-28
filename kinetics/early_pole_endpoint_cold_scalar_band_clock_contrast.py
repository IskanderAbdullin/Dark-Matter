"""Compare known same-restart S2/S4 q-band integrals at both final time steps."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_scalar_band_ladder import verify_pair


def contrasts(a,b,c,d):
    """a,b = S2 coarse/fine; c,d = S4 coarse/fine; all arrays use the same band."""
    a,b,c,d=[np.asarray(t,dtype=float) for t in (a,b,c,d)]
    if len({t.shape for t in (a,b,c,d)})!=1 or not np.all(np.isfinite([a,b,c,d])) or np.min([a,b,c,d])<=0:
        raise ArithmeticError('Need positive resolved coefficients on both clocks')
    coarse,fine=c/a-1,d/b-1
    return dict(S2_coarse_fine_minus_one=(a/b-1).tolist(),
                S4_coarse_fine_minus_one=(c/d-1).tolist(),
                S4_over_S2_on_coarse_clock_minus_one=coarse.tolist(),
                S4_over_S2_on_fine_clock_minus_one=fine.tolist(),
                change_in_contrast_fine_minus_coarse=(fine-coarse).tolist())


def run(files,states,archive):
    levels=(2,2,4,4)
    steps=(1024,2048,1024,2048)
    rows=[verify_pair(p,s,archive,l,steps=k) for p,s,l,k in zip(files,states,levels,steps,strict=True)]
    cols=[r['diagnostic'] for r in rows]
    vals=[c['raw_number_scalar_energy_and_excess'] for c in cols]
    eta=[c['conditional_number_energy_neutral'] for c in cols]
    poles=[c['full_pole_and_snapshot'] for c in cols]
    p=[[f[k] for k in ('scalar_mass','omega','residue')] for f in poles]
    return dict(grid='fine',x=.05,band=cols[0]['band'],charge_subcells=2,
                charts=levels,step_counts=steps,
                full_state_panel_receipts=rows,coefficient_contrasts=contrasts(*vals),
                neutral_root_contrasts=contrasts(*eta),optical_parameter_contrasts=contrasts(*p),
                scope='Paired time sensitivity for the final .04-to-.05 interval with each chart using its own fixed archived .04 state. Does not infer full-history time convergence or a physical source.')


if __name__=='__main__':
    r=argparse.ArgumentParser(description=__doc__)
    for name in ('s2-coarse','s2-fine','s4-coarse','s4-fine'):
        r.add_argument('--'+name,type=Path,nargs=2,required=True)
    r.add_argument('--states',type=Path,nargs=4,required=True)
    r.add_argument('--archive',type=Path,required=True)
    a=r.parse_args()
    a=(run((a.s2_coarse,a.s2_fine,a.s4_coarse,a.s4_fine),a.states,a.archive))
    print(json.dumps(a,indent=2,allow_nan=False))
