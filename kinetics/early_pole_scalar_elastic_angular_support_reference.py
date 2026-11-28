"""Kinematic angular support of elastic final states, not a Bose transport solver."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


def support(s,order):
    r=s['checkpoints'][-1]
    q=np.asarray(s['scalar_momentum_nodes'])
    n=np.maximum(np.asarray(s['endpoint_scalar_signed_shell_number_over_T3']),0.)
    edges=np.asarray(s['scalar_momentum_edges'])
    ir=edges[1:]<=1e-5
    f=np.asarray(s['endpoint_scalar_occupations'])
    pk=int(np.argmax(f))
    k0,k1=edges[pk:pk+2]
    m=r['scalar_pole_over_T']
    E=np.hypot(m,q)
    p=q[:,None]
    k=q[None,:]
    a=E[:,None]; b=E[None,:]
    qp=q*q/(E+m)
    C=(qp[:,None]+qp[None,:])/2
    d0=(p*p*b/(a+m)+m*k*k/(b+m))
    wt=n[:,None]*n[None,:]
    incoming=ir[:,None]&ir[None,:]
    total=0.
    cols=np.zeros(5)
    zz,ww=roots_legendre(order)
    KC=lambda z:z*z/(np.hypot(m,z)+m)
    for z,u in zip(zz,ww,strict=True):
        d=np.maximum(d0-p*k*z,0.)
        shat=4*m*m+2*d
        lab=np.sqrt(np.maximum(p*p+k*k+2*p*k*z,0.))
        vg=np.sqrt(d*(2*m*m+d))/(a*b)
        T=lab*np.sqrt(d/(2*shat))
        w=wt*vg/shat*u
        total+=w.sum()
        def both(K1,K2):
            gap=np.minimum(C-K1,K2-C)
            # If boost vanishes, both daughters have identical energies.
            raw=np.divide(gap,T,out=np.where(gap>=0.,1.,0.),where=(T>1e-30))
            return np.clip(raw,0.,1.)
        v=w*incoming
        cols+=np.array([v.sum(),np.sum(v*both(0.,KC(1e-5))),
                        np.sum(v*both(KC(k0),KC(k1))),
                        np.sum(v*both(0.,KC(.1))),
                        np.sum(w*both(0.,KC(1e-5)))])
    assert total>0 and cols[0]>0
    if np.max(q[ir])**2*2 < .1**2:
        assert abs(cols[3]/cols[0]-1)<1e-12
    return dict(q_peak_cell=[float(k0),float(k1)],
                fraction_incoming_both_q_below_1e_5=float(cols[0]/total),
                outgoing_both_q_below_1e_5_given_incoming_ir=float(cols[1]/cols[0]),
                outgoing_both_in_peak_cell_given_incoming_ir=float(cols[2]/cols[0]),
                outgoing_both_q_below_point1_given_incoming_ir=float(cols[3]/cols[0]),
                outgoing_both_q_below_1e_5_total=float(cols[4]/total))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spectra',nargs='+')
    args=p.parse_args()
    out={}
    for name in args.spectra:
        f=Path(name)
        raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        d=json.loads(raw)
        first,last=support(d,64),support(d,128)
        dr=max(abs(last[k]-first[k]) for k in last if k!='q_peak_cell')
        assert dr<2e-4,(f.name,dr)
        last['max_angular_refinement_drift']=dr
        out[f.name]=last
    print(json.dumps(out,indent=2,allow_nan=False))
