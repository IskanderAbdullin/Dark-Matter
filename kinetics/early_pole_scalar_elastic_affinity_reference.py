"""Stable 2->2 Bose collision affinity, including fractional energy stencils."""
import json

import numpy as np


def stable_gain(f,w):
    """Net forward occupation factor; w sums to two on each side.

    Positive finite occupations are required, including for fractional cells.
    This is a static algebra test, not a scalar transport solver.
    """
    f,w=np.asarray(f),np.asarray(w)
    assert f.shape==w.shape and f.shape[0]==2
    assert np.all(f>0.) and np.all(w>=0.) and np.allclose(w.sum(axis=1),2.)
    log=np.log(f)
    Bose=np.logaddexp(0.,-log)
    A=np.sum(w*Bose,axis=1)
    D=A[1]-A[0]
    scale=np.exp(np.sum(w*log)+min(A))
    return scale*np.sign(D)*np.expm1(abs(D))


def check():
    f=np.array([[2e16,4e16],[1e16,7e16]])
    w=np.ones_like(f)
    naive=f[0].prod()*np.prod(1+f[1])-f[1].prod()*np.prod(1+f[0])
    a,b=f[0]
    c,d=f[1]
    target=a*b*(1+c+d)-c*d*(1+a+b)
    value=stable_gain(f,w)
    assert naive==0 and abs(value/target-1)<1e-12
    assert abs(stable_gain(f[::-1],w[::-1])/value+1)<1e-12
    # Incoming particles each at a node; two outgoing particles split
    # across neighbors with common total average energy.
    E=np.array([[2.,5.,1.,1.],[2.,4.,3.,6.]])
    W=np.array([[1.,1.,0.,0.],[.5,.5,2/3,1/3]])
    assert abs(W[0]@E[0]-W[1]@E[1])<1e-14
    F=1./np.expm1(1e-17*E)
    q=stable_gain(F,W)
    scale=np.prod(F.ravel()**W.ravel())*1e-17
    assert abs(q/scale)<3e-13
    return dict(high_occupation_direct_rate=float(naive),high_occupation_stable_rate=float(value),
                signed_cubic_relative_error=float(abs(value/target-1)),
                normalized_fractional_thermal_balance_residual=float(abs(q/scale)))


if __name__=='__main__':
    print(json.dumps(check(),indent=2,allow_nan=False))
