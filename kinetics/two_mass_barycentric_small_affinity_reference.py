"""Small Bose/Pauli affinity with separate dominant and sub-ulp shell logits."""
import json

import numpy as np

from two_mass_saturated_bose_gradient_reference import DT, bose_average


def af(parent,fermion,bose_exp,wb,wf):
    p=np.asarray(parent,dtype=DT)
    q=np.asarray(fermion,dtype=DT)
    z=np.asarray(bose_exp,dtype=DT)
    w=np.asarray(wb,dtype=DT)
    v=np.asarray(wf,dtype=DT)
    if min(w)<=0 or min(v)<=0 or min(z)<=0 or q.shape!=(len(v),2) or z.shape!=w.shape:
        raise ValueError("Positive weights, Bose variables and two-part Fermi exponents required.")
    w,v=w/w.sum(),v/v.sum()
    a=w@z
    b=v@q[:,0]+v@q[:,1]
    c=p[0]+p[1]
    A=a+v@((q[:,0]-p[0])+(q[:,1]-p[1]))
    fwd=-np.logaddexp(DT(0),c)-np.logaddexp(DT(0),-b)-np.log(-np.expm1(-a))
    out=np.sign(A)*np.exp(max(fwd,fwd-A))*(-np.expm1(-abs(A)))
    return A,out,a+b-c


def mp_flux(parent,fermion,bexp,wb,wf):
    import mpmath as mp
    def conv(x):
        n,d=DT(x).as_integer_ratio()
        return mp.mpf(n)/d
    with mp.workdps(100):
        p=sum(map(conv,parent))
        x,w,v=[list(map(conv,q)) for q in (bexp,wb,wf)]
        a=sum(y*t for y,t in zip(x,w,strict=True))/sum(w)
        f=sum(y*sum(map(conv,t)) for y,t in zip(v,fermion,strict=True))/sum(v)
        B=1/mp.expm1(a)
        F=1/(1+mp.exp(f))
        P=1/(1+mp.exp(p))
        return float(P*(1-F)*(1+B)-(1-P)*F*B)


def run():
    # A tiny but explicitly represented Bose exponent on two receiving
    # shells, with an additional perturbation smaller than the ulp of
    # either absolute Fermi exponent.
    z=[bose_average(DT(f),DT(f)*DT(q),DT(1))[0] for f,q in ((2e22,1.004),(5e22,.995))]
    wb,wf=[DT(.32),DT(.68)],[DT(.57),DT(.43)]
    v=np.array(wf,dtype=DT); v/=sum(v)
    w=np.array(wb,dtype=DT); w/=sum(w)
    offset=np.array([DT(-6e-23),DT(9e-23)])
    f=[[DT(.35),x] for x in offset]
    out=[]
    for request in (DT(-1e-28),DT(+1e-28),DT(3e-29)):
        parent=[DT(.35),w@z+v@offset-request]
        A,S,raw=af(parent,f,z,wb,wf)
        ref=mp_flux(parent,f,z,wb,wf)
        err=float(abs(DT(ref)/S-1))
        if A*request<=0 or err>2e-8 or abs(A/request-1)>5e-10:
            raise ValueError("Split exponent test did not resolve the intended small signed flux.")
        if raw!=0.:
            raise ValueError("The manufactured direct absolute exponent did not lose the tiny perturbation.")
        out.append(dict(requested_affinity=float(request),calculated_affinity=float(A),output_flux=float(S),
                        unsplit_direct_affinity=float(raw),mp_independent_flux=ref,relative_flux_error=err))
    return dict(scope="Manufactured high-Bose arithmetic, no full nonlinear saturated collision history.",cases=out)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
