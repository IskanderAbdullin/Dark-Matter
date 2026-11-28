"""Check fitted-rate replay and finite witness signs in the occupation envelope."""
import argparse
import json
from pathlib import Path


def compare(reference,runs):
    if len(reference)!=4 or len(runs)!=4 or any(len(v)!=4 for v in runs):
        raise ValueError('All four frozen aperture fixtures and constraint controls are required.')
    out=[]
    for a,*rows in zip(reference,*runs,strict=True):
        depth=a['hard_depth']
        modes=a['cold_modes']
        check=[]
        for v,b in zip(rows,[(False,4),(False,12),(True,4),(True,12)],strict=True):
            if (v['hard_depth'],v['cold_modes'],v['fixed_gauge_moment'],v['subcells_per_packet_interval'])!=(depth,modes,*b):
                raise ValueError('Mismatched fixed support or constraint set.')
            diff=abs(v['independent_fitted_rate']-a['reconstructed_paired_total'])
            internal=v['piecewise_witness_range']
            outer=v['integrated_dual_range']
            if diff>1e-11 or not outer[0]<=internal[0]<=internal[1]<=outer[1]:
                raise ValueError('Independent rate or extremal occupation check failed.')
            percell=[c[k] for c in v['cells'] for k in ('native','archive_local')]
            check.append(dict(fixed_gauge=v['fixed_gauge_moment'],refinement=v['subcells_per_packet_interval'],
                              fitted_reference_difference=diff,witness_range=internal,dual_range=outer,
                              spans_both_signs=internal[0]<0.<internal[1],
                              max_abs_moment_residual=max(c['max_abs_moment_residual'] for c in percell),
                              max_dual_quadrature_order_delta=max(c['max_dual_order_delta'] for c in percell)))
        out.append(dict(hard_depth=depth,cold_modes=modes,controls=check))
    return dict(scope='Explicit positive piecewise constant witnesses and exploratory integrated continuous dual evaluations; no propagated dynamics.',
                checked=out,raw=runs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference')
    p.add_argument('receipts',nargs=4)
    v=p.parse_args()
    rd=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(rd(v.reference),[rd(p) for p in v.receipts]),indent=2,allow_nan=False))
