"""Look for paired numerical counterexamples to interpreting a split as an error bound."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from early_pole_endpoint_scalar_subshell_spectrum_audit import read
from early_pole_endpoint_stencil_remap_pilot import load
from early_pole_endpoint_warm_refinement_transfer_test import source


def limits(q):
    delta=abs(q['split']-q['base'])/q['scale']
    one=np.asarray([q['reference32'],q['reference']])
    err=abs(one-q['base'])/q['scale']
    new=abs(one-q['split'])/q['scale']
    missed=(err-delta).min(axis=0)
    worse=(new-err).min(axis=0)
    def result(scores):
        j=int(np.argmax(scores))
        return dict(edges=q['edges'][j:j+2].tolist(),excess_under_both_references=float(scores[j]),
                    estimate=float(delta[j]),unrefined_L1_contribution=err[:,j].tolist(),
                    refined_L1_contribution=new[:,j].tolist())
    return dict(x=q['x'],rank_correlation=[float(spearmanr(delta,u).statistic) for u in err],
                number_underestimated_by_both=int(np.count_nonzero(missed>0)),
                number_worse_against_both_after_split=int(np.count_nonzero(worse>0)),
                strongest_underestimate=result(missed),strongest_degradation_after_split=result(worse))


def run(cards,refs):
    if len(cards)!=4 or [x['x_end'] for x in cards]!=[.01,.02,.04,.05]:
        raise ArithmeticError('Need the four matching archived warm snapshots')
    c=load('fine')[0]
    return dict(rows=[limits(source(x,c,r)) for x,r in zip(cards,refs,strict=True)],reference_order=[32,64],
                scope='Both references finite: these are numerical counterexamples to a local surrogate bound, not rigorous continuum error intervals.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--states',type=Path,nargs=4,required=True)
    p.add_argument('--references',type=Path,nargs=4,required=True)
    a=p.parse_args()
    print(json.dumps(run([read(t) for t in a.states],[read(t) for t in a.references]),indent=2,allow_nan=False))
