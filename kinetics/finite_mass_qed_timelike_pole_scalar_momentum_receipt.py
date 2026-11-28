"""Verify the same photon cannot balance the two archived scalar modes."""
import argparse
import json
from pathlib import Path


def run(paths):
    out={}
    bound=[]
    for g,path in zip(('fine','third'),paths,strict=True):
        vals=[json.loads(z) for z in path.read_text().splitlines()]
        v=[x for x in vals if x['pol']=='T' and x['k']==.0015 and x['q']<=.0008]
        expected=(0.,2.05e-6,4.1e-6,8.2e-6,4e-5,1e-4,3e-4,.0008)
        if tuple(x['q'] for x in v)!=expected or any(x['grid']!=g or x['mesh_defect']>3e-6 for x in v):
            raise ArithmeticError(('Mismatched scalar momentum ladder',g))
        q0,q1=v[2],v[-1]
        low,high=q0['isotropic_conditional_population'],q1['isotropic_conditional_population']
        n=1e-4  # comparison point, not a physical population estimate
        if not 0<low<n<high:
            raise ArithmeticError('Nonempty conditional interval not verified')
        J=[]
        for x in (q0,q1):
            A,B,D=x['raw_A_B_D']
            J.append(n*(D+A/x['scalar_occup'])-B)
        if not J[0]>0>J[1]:
            raise ArithmeticError('Signed affinity type not verified')
        bound.extend([(low,high)])
        out[g]=dict(samples=v,trial_population=n,physical_q_pair=[q0['q'],q1['q']],
                    opposite_signs_before_coupling=[float(x) for x in J],
                    threshold_ratio=high/low)
    return dict(scope='Conditional fixed transverse photon k=.0015 only. No time-evolved pole population.',
                shared_opposite_sign_open_interval=[max(x[0] for x in bound),min(x[1] for x in bound)],
                endpoints=out)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('files',type=Path,nargs=2)
    print(json.dumps(run(p.parse_args().files),indent=2,allow_nan=False))
