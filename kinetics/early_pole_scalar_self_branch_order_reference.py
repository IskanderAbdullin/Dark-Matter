"""Reversing the order of matched screened comparisons must not change results."""
import hashlib
import json

from early_pole_joined_scalar_self_proxy_reference import run


def result(reverse):
    a=run('fine',16,6,2e-9,2,3,hard_slow=True,screened=True,mixes=(.97,),reverse_branches=reverse)
    return sorted(a['cases'],key=lambda g:g['scalar_elastic_proxy'])


if __name__=='__main__':
    x,y=result(False),result(True)
    h=lambda z:hashlib.sha256(json.dumps(z,allow_nan=False,sort_keys=True).encode()).hexdigest()
    assert x==y
    print(json.dumps(dict(forward_hash=h(x),reverse_hash=h(y),cases=[
        dict(scalar_proxy=r['scalar_elastic_proxy'],gain=r['rows'][-1]['scalar_gain']) for r in x]),
        indent=2,allow_nan=False))
