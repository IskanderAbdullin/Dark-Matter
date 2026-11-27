"""Keep both refinement arms visible when their paired difference is small."""
import argparse
import gzip,json
from math import fsum
from pathlib import Path

import numpy as np


def read(p):
    if p.endswith('.gz'):
        with gzip.open(p,'rt') as f:return json.load(f)
    return json.loads(Path(p).read_text())


def analyze(low,high,low_depth,high_depth,check=None):
    z=[[read(f) for f in pair] for pair in (low,high)]
    v=z[0][0]
    coords=['projection','rtol','base_atol','duration','history_exposure','cold_momentum_nodes_over_T','report_panels']
    T=np.asarray(v['history_exposure'])
    a=np.asarray(v['cold_mode_count_history'])[0]
    history=[]
    for pair,d in zip(z,(low_depth,high_depth),strict=True):
        buf=[]
        for arm,r in zip(('native','local-daughter'),pair,strict=True):
            observed=(r.get('aperture_subdivision',1),r.get('prepared_charged_swap','native'),
                      r.get('selected_elastic_channel','both'))
            if (observed!=(d,arm,'both') or not r.get('elastic_channels',True)
                    or any(r[k]!=v[k] for k in coords) or abs(r['energy_defect'])>=3e-12
                    or abs(r['entropy_defect'])>=3e-11 or abs(r['hard_number_defect'])>=3e-12):
                raise ArithmeticError('The arms did not come from the matched held reference.')
            n=np.asarray(r['cold_mode_count_history'])
            if n.shape!=(len(T),32) or not np.all(np.isfinite(n)&(n>0)):
                raise ArithmeticError('Incomplete cold modes.')
            np.testing.assert_array_equal(n[0],a)
            np.testing.assert_array_equal(n.T.copy(order='C').sum(axis=0),r['cold_count_history'])
            buf.append(n-a)
        history.append(buf)
    dx=np.diff(np.asarray(history),axis=0)[0] # arm, time, mode
    δ=np.array([[fsum(y) for y in x] for x in dx])
    reshape=np.diff(dx,axis=1)
    dT=np.diff(T)
    rate=np.array([[fsum(y)/dt for y,dt in zip(x,dT,strict=True)] for x in reshape])
    paired=np.array([fsum(x) for x in dx[1]-dx[0]])
    closes=max(abs(δ[1]-δ[0]-paired))
    if closes>2e-22:
        raise ArithmeticError('Direct arm contrasts did not recover the modal paired correction.')
    if check:
        np.testing.assert_allclose(paired,check['cumulative_h8_minus_h4'],rtol=0.,atol=2e-22)
        np.testing.assert_allclose(rate[1]-rate[0],check['adjacent_panel_mean_rates'],rtol=0.,atol=1e-14)
    common=rate.mean(axis=0)
    Paired=rate[1]-rate[0]
    report=[]
    for i,(n,l) in enumerate(zip(*rate,strict=True)):
        amount=abs(n)+abs(l)
        report.append(dict(panel=i+1,midpoint=float((T[i]+T[i+1])/2),
                           native_refinement_rate=float(n),local_refinement_rate=float(l),
                           local_minus_native=float(Paired[i]),
                           mean_common_refinement_rate=float(common[i]),
                           both_individual_rate_signs_negative=bool(n<0 and l<0),
                           native_positive_modes=int(np.count_nonzero(reshape[0,i]>0)),
                           local_positive_modes=int(np.count_nonzero(reshape[1,i]>0)),
                           signed_cancellation_fraction=float(1-abs(Paired[i])/amount) if amount else None,
                           common_to_paired_magnitude_ratio=float(abs(common[i]/Paired[i])) if Paired[i] else None))
    s=abs(δ[:,-1]).sum()
    return dict(scope='Separate H4→H8 arms inside the held finite archived control. No extrapolation.',
                depth_change=[low_depth,high_depth],exposure=T.tolist(),native_and_local_cumulative=δ.tolist(),
                paired_cumulative=paired.tolist(),terminal_arm_differences=δ[:,-1].tolist(),
                terminal_negative_modes=np.count_nonzero(dx[:,-1]<0,axis=-1).tolist(),
                terminal_cancellation_fraction=float(1-abs(paired[-1])/s) if s else None,panel_readout=report)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--low',nargs=2,required=True)
    p.add_argument('--high',nargs=2,required=True)
    p.add_argument('--low-depth',type=int,default=4)
    p.add_argument('--high-depth',type=int,default=8)
    p.add_argument('--paired-receipt')
    a=p.parse_args()
    print(json.dumps(analyze(a.low,a.high,a.low_depth,a.high_depth,read(a.paired_receipt) if a.paired_receipt else None),
                     indent=2,allow_nan=False))
