"""Stable conditional-threshold decomposition of previously measured frozen scalar points."""
import argparse
import json
from pathlib import Path


def run(data):
    c0,c1=data['cards']
    if (c0['scalar_div'],c1['scalar_div'])!=(2,4):
        raise ArithmeticError('Expected recorded S2→S4 fixed-q points')
    out=[]
    for l,r in zip(c0['points'],c1['points'],strict=True):
        if l['q']!=r['q']:
            raise ArithmeticError('Incompatible q')
        A0,B0,D0=l['raw_affine']
        A1,B1,D1=r['raw_affine']
        f0,f1=l['occupation'],r['occupation']
        if min(A0,B0,D0,A1,B1,D1,f0,f1)<=0:
            raise ArithmeticError('The stable positive-balance expression does not apply')
        # Stabilized form, even if rounding would make n(f1)/n(f0) equal one.
        f_only=(A0/D0)*(1/f0-1/f1)/(1+A0/(D0*f1))
        kernel_only=((B1-B0)*(D0+A0/f1)-B0*((D1-D0)+(A1-A0)/f1))/(B0*(D1+A1/f1))
        actual=r['conditional_photon_number_balance']/l['conditional_photon_number_balance']-1
        if abs(f_only+kernel_only+f_only*kernel_only-actual)>2e-14:
            raise ArithmeticError('Sequential balance accounting did not close')
        out.append(dict(q=l['q'],S2_scalar_only_fractional_shift=f_only,
                        subsequent_raw_kernel_fractional_shift=kernel_only,combined_fractional_shift=actual))
    return dict(cases=out,scope='Fixed-q positive conditional crossings; actual occupation not determined.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt',type=Path)
    print(json.dumps(run(json.loads(p.parse_args().receipt.read_text())),indent=2,allow_nan=False))
