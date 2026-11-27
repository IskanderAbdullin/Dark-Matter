"""Finite triangular eikonal current, exact Abelian conservation check."""
import json

import numpy as np


def segment(K,x,y):
    d=y-x
    angle=K@d
    return np.exp(1j*(K@((x+y)/2)))[:,None]*np.sinc(angle/(2*np.pi))[:,None]*d


def run():
    rng=np.random.default_rng(427)
    K=np.r_[rng.normal(size=(4096,4)),np.eye(4),np.zeros((1,4))]
    out=[]
    for t in (.02,1.,30.,240.):
        for dv in (0.,.0006036,.03):
            root=np.zeros(4)
            P=t*np.array([1.,0.,0.,.9983])
            B=t*np.array([1.,0.,0.,.9983+dv])
            p,b,c=segment(K,root,P),segment(K,B,root),segment(K,P,B)
            closed=(p+b)+c
            a=np.einsum('ij,ij->i',K,closed)
            op=np.einsum('ij,ij->i',K,p+b)
            expectation=(np.exp(1j*K@P)-np.exp(1j*K@B))/1j
            norm=max(float(np.linalg.norm(closed,axis=1).max()),1.)
            residual=max(abs(a))/norm
            assert residual<3e-12 and max(abs(op-expectation))<3e-12
            if dv==0:
                assert np.max(np.abs(closed))==0.
            out.append(dict(t=t,dv=dv,expanded_current_norm_max=norm,
                            closed_cancellation_relative_max=residual,
                            lost_connector_open_contraction_max=max(abs(op))))
    return dict(scope='External prescribed paths only; no fitted noise, rates, bath, or equilibration.',rows=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
