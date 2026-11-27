"""Describe successive first-window refinements without claiming a continuum limit."""
import argparse
import json
import math

from early_pole_joined_initial_window_refinement_reference import scores
from early_pole_joined_temporal_refinement_reference import select


def cauchy(paths,projection):
    runs=[select(p,projection) for p in paths]
    result=scores(runs)
    div=result['resolution_divisors']
    if len(div)<3 or not div[-1]==2*div[-2]==4*div[-3]:
        raise ValueError("The last three estimates must be on dyadic widths.")
    a,b,c=result['first_interval_averages'][-3:]
    d1,d2=a-b,b-c
    observed=math.log2(abs(d1/d2)) if d1*d2>0 else None
    result.update(comparison_dyadic_divisors=div[-3:], signed_differences=[d1,d2],
                  observed_local_order=observed,
                  linear_extrapolation_diagnostic=c-d2 if observed is not None else None)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='+')
    p.add_argument('--projection',choices=('steps','cell-conservative'),default='steps')
    a=p.parse_args()
    print(json.dumps(cauchy(a.paths,a.projection),indent=2,allow_nan=False))
