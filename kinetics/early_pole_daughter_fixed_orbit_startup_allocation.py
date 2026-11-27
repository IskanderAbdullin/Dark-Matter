"""Localize the fixed-inner time gap before evolving a slow-coupled hybrid."""
import argparse
import json
from itertools import combinations
from math import fsum
from pathlib import Path

import numpy as np


def compare(orbit, dyadic, hybrid1, hybrid2):
    q=np.asarray(orbit['cold_source_by_mode'])
    if q.shape!=(2,65,8) or not np.all(np.isfinite(q)):
        raise ValueError('Require both held orbits and 65 complete cold-source samples.')
    full=dyadic['results']['16']['individual_outer_comparison']
    if orbit['leaf_width']!=full['exposure'][-1]/64:
        raise ValueError('The elastic leaves differ.')
    nodes=np.arange(16)[:,None]*4+np.arange(5)
    win=q[:,nodes,:]
    # Each block contributes one-sixteenth of the final mean.
    native,local=np.einsum('abiq,i->abq',win,np.array([-3,2,2,2,-3])/128)
    pair=local-native
    per_block=pair.sum(axis=-1)
    total=fsum(per_block)
    static=orbit['cases']
    t0=static['16']['original_trapezoid_mean']
    tf=static['32']['refined_trapezoid_mean']
    if not np.isclose(total,tf-t0,rtol=0.,atol=2e-17):
        raise ValueError('Selected block contributions do not close to the fixed-orbit total.')
    cases=full['cases']
    coarse=cases['elastic_centered']['final_mean_shift']
    fine=cases['fine64']['final_mean_shift']
    signed_gap=fine-coarse

    ranked={}
    for k,v in (('native_total',np.abs(native.sum(axis=-1))),
                ('local_total',np.abs(local.sum(axis=-1))),
                ('paired_total',np.abs(per_block)),
                ('sum_of_all_arm_mode_magnitudes',np.abs(local).sum(axis=-1)+np.abs(native).sum(axis=-1))):
        ranked[k]=(np.argsort(-v,kind='stable')+1).tolist()
    optimum={}
    for n in range(1,5):
        selected=min(combinations(range(16),n),
                     key=lambda pick:abs(total-fsum(per_block[j] for j in pick)))
        optimum[str(n)]=dict(blocks=[j+1 for j in selected],
                             residual_mean=float(total-fsum(per_block[j] for j in selected)))
    data={}
    for n,other in ((1,hybrid1),(2,hybrid2)):
        table=other['cases']
        if (table['coarse16x4']['mean_shift']!=coarse or table['fine64x1']['mean_shift']!=fine
                or table[f'hybrid{16+3*n}']['outer_intervals']!=16+3*n):
            raise ValueError('Historical hybrid measurements use a different comparison.')
        actual=table[f'hybrid{16+3*n}']['mean_shift']
        shift=fsum(per_block[:n])
        data[str(n)]=dict(static_recovered_fraction=shift/total,observed_recovered_fraction=(actual-coarse)/signed_gap,
                          observed_mean=actual,static_plus_observed_coarse_mean=coarse+shift,
                          observed_minus_prediction=actual-coarse-shift)
    later=float(per_block[1])
    return dict(scope='Offline block scores use all 65 fixed-orbit samples; they are not an online error bound or a new run of a slow-coupled trajectory.',
                block_number=list(range(1,17)),
                native_mean_correction_by_mode=native.tolist(), local_mean_correction_by_mode=local.tolist(),
                paired_mean_correction_by_mode=pair.tolist(), paired_mean_correction_by_block=per_block.tolist(),
                rankings=ranked, optimal_small_budget_for_this_paired_orbit=optimum,
                first_to_first16_static_recovered_fraction=(np.cumsum(per_block)/total).tolist(),
                fixed_orbit_full_gap=total, slow_coupled_full_gap=signed_gap,
                previous_hybrid_controls=data,
                second_only_predeclared_prediction=dict(paired_mean=coarse+later,
                                                        corrected_fraction_of_static_gap=later/total,
                                                        refined_blocks=[2]))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('orbit','dyadic','hybrid1','hybrid2'):
        p.add_argument(k)
    v=p.parse_args()
    print(json.dumps(compare(*(json.loads(Path(f).read_text()) for f in vars(v).values())),
                     indent=2,allow_nan=False))
