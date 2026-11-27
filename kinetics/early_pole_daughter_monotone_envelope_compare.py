"""Cross-check monotone cutoff envelopes against the original fitted source."""
import argparse
import json
from pathlib import Path


def compare(fitted,raw):
    if len(fitted)!=4 or len(raw)!=4 or any(len(v)!=4 for v in raw):
        raise ValueError('Need four accepted hard/cold fixture controls.')
    out=[]
    for f,*allv in zip(fitted,*raw,strict=True):
        key=(f['hard_depth'],f['cold_modes'])
        row=[]
        for a,kind in zip(allv,[(False,12),(False,48),(True,12),(True,48)],strict=True):
            if (a['hard_depth'],a['cold_modes'])!=key or (a['fixed_gauge_moment'],a['subdivision'])!=kind:
                raise ValueError('Unexpected selected fixture or cutoff refinement.')
            actual=a['independent_fitted_pair']
            err=abs(actual-f['reconstructed_paired_total'])
            x,y=a['numerical_outer_paired_interval'],a['relaxed_native_numerical_outer_interval']
            if err>1e-11 or not y[0]<=x[0]<=actual<=x[1]<=y[1]:
                raise ValueError('Pre-existing fitted control or relaxed envelope did not match.')
            d=[v[k]['dual_order_32_to_64_delta'] for v in a['cells'] for k in (
                'local_cutoff_min','local_cutoff_max','native_relaxed_min','native_relaxed_max')]
            row.append(dict(fixed_gauge_moment=kind[0],subdivision=kind[1],
                            outer_interval=x,relaxed_native_outer_interval=y,
                            max_abs_dual_quadrature_delta=max(d),fitted_reference_error=err,
                            max_abs_native_flat_moment_defect=max(v['flat_native_energy_moment_defect'] for v in a['cells'])))
        out.append(dict(hard_depth=key[0],cold_modes=key[1],controls=row))
    return dict(scope='Conditional monotonicity with frozen packet geometry and numerical cutoff dual integrals.',
                matched_fixtures=out,raw=raw)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference')
    p.add_argument('receipts',nargs=4)
    a=p.parse_args()
    load=lambda f:json.loads(Path(f).read_text())
    print(json.dumps(compare(load(a.reference),[load(f) for f in a.receipts]),indent=2,allow_nan=False))
