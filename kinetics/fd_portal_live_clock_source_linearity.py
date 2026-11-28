"""Check local linearized source error at independently evolved FD clock states."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from fd_portal_combined_tabulated_bath_entropy import elastic_heat
from fd_portal_differential_reciprocity_reference import PositivePortal,integrate


def check(p,mesh,T,eta,beta,combined):
    g=np.array([2*eta,beta])
    s=integrate(p,T,*g)[0]
    _,J,_,_=integrate(p,T,0.,0.)
    a=s[:2]*[2.,T]
    b=(J@g)*[2.,T]
    if combined:
        a[1]+=elastic_heat(T,T/(1-beta),eta,mesh)
        b[1]+=-elastic_heat(T,T,0.,mesh,linear=True)*beta
    m=np.array([eta,beta/T])
    if (m@a)<=0. or (m@b)<=0. or min(abs(a))<=1e-90:
        raise ArithmeticError('Resolved physical and linearized dissipative sources required')
    vals=np.r_[a,m@a]
    linear=np.r_[b,m@b]
    return vals,linear


def run(archive):
    z=archive.read_bytes(); d=json.loads(z)
    p,q=[PositivePortal(o) for o in ((13,36,26),(17,44,32))]
    rows=[]
    for r in d['rows']:
        T,eta,beta=[r[k] for k in ('local_crossing_T','actual_eta','actual_beta')]
        u,v=[check(i,n,T,eta,beta,r['operator']=='combined')
             for i,n in ((p,(18,29.)),(q,(25,36.)))]
        err=np.max(abs(np.array(u)/v-1.))
        if err>3e-5:
            raise ArithmeticError('Clock-specific actual source did not resolve under both phase-space rules')
        a,b=v
        rows.append(dict(operator=r['operator'],mode_label=r['mode'],T_SM=T,eta_and_beta=[eta,beta],
                         direct_number_heat_entropy=a.tolist(),linear_number_heat_entropy=b.tolist(),
                         linear_over_direct_minus_one=(b/a-1.).tolist(),same_sign=[bool(x*y>0.) for x,y in zip(a,b,strict=True)],
                         maximum_relative_mesh_difference=err))
    if len(rows)!=4:
        raise ArithmeticError('Require both numerical clocks for both direct histories')
    return dict(rows=rows,archive_sha256=hashlib.sha256(z).hexdigest(),
                scope='Instantaneous source comparison on resolved finite states, not a stationary or relic-abundance model.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    print(json.dumps(run(p.parse_args().archive),indent=2,allow_nan=False))
