"""Bound what a measured cold response can imply once an arbitrary unmeasured tail is added."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np


def tv_interval(p,q,e,f):
    p,q=np.asarray(p,dtype=float),np.asarray(q,dtype=float)
    if (p.shape!=q.shape or p.ndim!=1 or min(p)<0 or min(q)<0 or
        max(abs(fsum(p)-1),abs(fsum(q)-1))>4e-14 or not 0<=e<=1 or not 0<=f<=1):
        raise ArithmeticError('Need two positive normalized cold profiles and valid missing response fractions')
    cold=fsum(abs((1-e)*p-(1-f)*q))
    return [(cold+abs(e-f))/2,(cold+e+f)/2]


def break_order(d_prev,d_next):
    if not 0<=d_next<d_prev<=1:
        raise ArithmeticError('Need a smaller later finite-band distance')
    change=d_prev-d_next
    return change/(1+change)


def worst_rank_gap(p,q,r,cap):
    """Maximum next TV minus previous TV if each omitted fraction is at most cap.

    At fixed middle fraction f, the smallest possible first distance has its
    left fraction equal to f. The next upper bound is nondecreasing in both
    of its fractions. Thus the simultaneous extreme takes all fractions=cap,
    with the first two exterior shapes identical and the third disjoint.
    """
    first=tv_interval(p,q,cap,cap)[0]
    second=tv_interval(q,r,cap,cap)[1]
    return second-first


def run(card):
    if card['levels']!=[2,4,8,16] or len(card['common_q_edges'])!=18:
        raise ArithmeticError('Use the published common 17-bin raw spatial profiles')
    results={}
    for field in ('forward_number','inverse_number','forward_excess'):
        profiles=[]
        totals=[]
        known=[]
        for z in card['adjacent']:
            o=z['observables'][field]
            x=np.asarray(o['native_and_crossed_bin_response_bounds'])
            if not np.array_equal(x[[0,1],[0,1],:,0],x[[0,1],[0,1],:,1]):
                raise ArithmeticError('Need exact native diagonal responses')
            rate=x[[0,1],[0,1],:,0]
            t=[fsum(v) for v in rate]
            val=rate/np.asarray(t)[:,None]
            d=tv_interval(*val,0.,0.)[0]
            if abs(d-o['native_profile_TV'])>2e-14:
                raise ArithmeticError('Raw native cold profile did not reproduce the archived distance')
            if profiles and (not np.array_equal(val[0],profiles[-1]) or t[0]!=totals[-1]):
                raise ArithmeticError('The middle chart changed between comparisons')
            if not profiles:
                profiles.append(val[0])
                totals.append(t[0])
            profiles.append(val[1])
            totals.append(t[1])
            known.append(d)
        rows=[]
        for i in (0,1):
            x,y=known[i:i+2]
            e=break_order(x,y)
            lower=tv_interval(profiles[i],profiles[i+1],e,e)[0]
            upper=tv_interval(profiles[i+1],profiles[i+2],e,e)[1]
            rank=worst_rank_gap(*profiles[i:i+3],e)
            if max(abs(upper-lower),abs(rank))>1e-14:
                raise ArithmeticError('Disjoint-versus-identical exterior witnesses did not meet')
            rows.append(dict(spatial_levels=card['levels'][i:i+3],
                             cold_native_tv_decrease=x-y,
                             minimum_common_outside_fraction_to_erase_contraction=e,
                             minimum_per_chart_omitted_fraction_cap_for_possible_reversal=e,
                             matched_distance_with_worst_outside_shapes=lower))
        amp=[]
        for i in range(3):
            a,b=totals[i:i+2]
            if not a>b>0:
                raise ArithmeticError('The retained integrated response is not a decreasing sequence')
            amp.append(dict(levels=card['levels'][i:i+2],cold_band_totals=[a,b],
                            extra_outside_response_finer_minus_coarser_to_erase_order=a-b,
                            extra_outside_response_relative_to_coarser_cold=(a-b)/a))
        results[field]=dict(cold_native_tv=known,shape_order_limits=rows,absolute_band_order_limits=amp)
    return dict(levels=card['levels'],common_q_edges=card['common_q_edges'],results=results,
                scope='Worst-case algebra for arbitrary exterior response. Analytical cap on each unknown omitted fraction, not measured tail masses, physical response or conditional profile optimum.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    a=p.parse_args()
    print(json.dumps(run(json.loads(a.source.read_text())),indent=2,allow_nan=False))
