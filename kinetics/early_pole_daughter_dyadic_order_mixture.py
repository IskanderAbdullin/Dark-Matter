"""Compare mixtures of opposite held orders and the next dyadic original order."""
import argparse
import json
from pathlib import Path

import numpy as np

from early_pole_daughter_alternate_outer_order_compare import compare


def validate(original16, alternate16, original32, alternate32, fine64):
    arms = {16:(original16, alternate16), 32:(original32, alternate32)}
    comparisons = {
        n:compare(dict(elastic_centered=pair[0], slow_centered=pair[1], fine64=fine64), n)
        for n,pair in arms.items()
    }
    seeds = [np.spacing(v['state_history'][0]) for v in fine64]
    refined_reference = {16:(32,'elastic_centered'), 32:(32,'fine64')}
    receipts = {}
    for n, result in comparisons.items():
        j,k = refined_reference[n]
        cases = result['cases']
        x,y = (np.asarray(cases[z]['paired_count_history']) for z in ('elastic_centered','slow_centered'))
        f = np.asarray(comparisons[j]['cases'][k]['paired_count_history'])
        error = (x+y)/2-f
        separation = x-y
        if not np.all(np.isfinite(error)):
            raise ValueError('The paired history contains non-finite values.')
        peak = float(max(abs(error)))
        disag = float(max(abs(separation)))
        if disag == 0:
            raise ValueError('Unresolved order contrast.')
        receipts[str(n)] = dict(
            individual_outer_comparison=result,
            refined_original_macro_steps=n*2,
            refined_original_paired_counts=f.tolist(),
            order_mixture_minus_refined_count=error.tolist(),
            peak_count_difference=peak, peak_coarse_order_separation=disag,
            peak_residual_over_peak_order_separation=peak/disag,
            peak_residual_over_largest_seed_spacing=peak/float(max(seeds)),
            final_mean_mixture_minus_refined=float(error[-1]/result['exposure'][-1]))
    return dict(scope='Two held outer orders at unchanged 64 elastic leaves; comparisons on the same 17 stored points. Arithmetic mean is not a single admissible evolved state.',
                results=receipts)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for tag in ('original16','alternate16','original32','alternate32','fine64'):
        p.add_argument('--'+tag, nargs=2, required=True, help='Native then local.')
    a = p.parse_args()
    print(json.dumps(validate(**{k:[json.loads(Path(v).read_text()) for v in vs]
                                 for k,vs in vars(a).items()}), indent=2, allow_nan=False))
