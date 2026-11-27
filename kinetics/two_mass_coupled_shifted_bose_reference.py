"""Manufactured coupled event limit with sub-ulp fermion count shifts."""
import json

import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.linalg import expm

DT=np.longdouble


def trial(raw_eps,f0,z0,c,parents,daughters,bose,time=.75):
    eps=DT(raw_eps)
    f=DT(f0)
    u=f*(1-f)
    z,c=np.asarray(z0,dtype=DT),np.asarray(c,dtype=DT)
    P,F,W=[np.asarray(v,dtype=DT) for v in (parents,daughters,bose)]
    for basis in (P,F,W):
        basis/=basis.sum(axis=0)
    V=np.r_[-P,F]
    theta=W.T@(1/c)
    if min(theta)<=0:
        raise ValueError("All event capacities and thermal exponents must be positive.")
    N=V.shape[1]

    def state(x):
        x=np.asarray(x,dtype=DT)
        vf=z+V@x
        vb=W@x
        if min(f+eps*vf)<=0 or max(f+eps*vf)>=1 or min(c+eps*eps*vb)<=0:
            raise ValueError("Episode left the positive Bose or Pauli interior.")
        delta=np.log1p(-eps*vf/(1-f))-np.log1p(eps*vf/f)
        pp=P.T@delta[:P.shape[0]]
        ff=F.T@delta[P.shape[0]:]
        alphaB=W.T@np.log1p(eps/(c+eps*eps*vb))
        fwdP=1/(1+(1-f)/f*np.exp(pp))
        noF=1/(1+f/(1-f)*np.exp(-ff))
        A=alphaB+ff-pp
        scaled_fwd=fwdP*noF*eps/(-np.expm1(-alphaB))
        S=scaled_fwd*(-np.expm1(-A))/eps
        return S,A/eps

    def vector(_,x):
        return np.asarray(state(x)[0],float)
    s0=np.zeros(N)
    sol=solve_ivp(vector,(0.,time),s0,method="DOP853",rtol=5e-13,atol=2e-15,dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    final=sol.y[:,-1]
    state_shift=z+V@final
    limit=np.zeros((N+1,N+1))
    limit[:N,:N]=np.asarray(-(V.T@V)/theta[:,None],float)
    limit[:N,-1]=np.asarray(u-(V.T@z)/theta,float)
    expected=expm(time*limit)[:N,-1]
    entropy_asym=theta@expected-(np.sum((z+V@expected)**2)-z@z)/(2*u)

    def rate(t):
        S,A=state(sol.sol(t))
        return float(S@A)
    a=quad(rate,0.,time,epsrel=2e-12,epsabs=3e-14)[0]

    import mpmath as mp
    def cv(x):
        n,d=DT(x).as_integer_ratio()
        return mp.mpf(n)/d
    with mp.workdps(135):
        ee,ff,zz,xx,cc=cv(eps),cv(f),list(map(cv,z)),list(map(cv,final)),list(map(cv,c))
        ss=lambda p:-p*mp.log(p)-(1-p)*mp.log1p(-p)
        sb=lambda n:mp.log1p(n)+n*mp.log1p(1/n)
        # Normalize each column *in the reference precision*. A residual
        # O(ulp) violation of exact number conservation would otherwise
        # pollute the entropy after dividing a linear term by eps².
        def mp_weight(block,sign):
            cols=list(zip(*block,strict=True))
            ps=[[cv(v) for v in q] for q in cols]
            cs=[[sign*v/sum(q) for v in q] for q in ps]
            return [sum(v*w for v,w in zip(q,xx,strict=True)) for q in zip(*cs,strict=True)]
        event_F=mp_weight(parents,-1)+mp_weight(daughters,1)
        event_B=mp_weight(bose,1)
        a0=[ff+ee*t for t in zz]
        b0=[t/ee for t in cc]
        mp_entropy=(sum(ss(p+ee*r)-ss(p) for p,r in zip(a0,event_F,strict=True))
                    +sum(sb(p+ee*r)-sb(p) for p,r in zip(b0,event_B,strict=True)))/ee**2

    match=abs(a/float(mp_entropy)-1)
    err=float(np.max(abs(final-expected)))
    nulls=[float(abs(np.sum(state_shift-z))),float(abs(sum((W@final))-sum(final)))]
    if match>4e-11 or a<=0. or max(nulls)>3e-15:
        raise ValueError(f"Independent entropy or number constraints failed: {raw_eps}, {match}, {nulls}.")
    if eps<=DT(1e-24) and (err>5e-13 or abs(a/entropy_asym-1)>5e-12):
        raise ValueError("Coupled high-occupancy limit did not match analytic affine flow.")
    ff=np.asarray(f+eps*z,float)
    ff2=np.asarray(f+eps*state_shift,float)
    return dict(epsilon=raw_eps,scaled_event_extents=final.tolist(),
                scaled_initial_signed_rates=vector(0,s0).tolist(),
                analytic_limit_extents=expected.tolist(),max_absolute_extent_error=err,
                integrated_scaled_entropy=a,analytic_limit_scaled_entropy=float(entropy_asym),
                high_precision_scaled_entropy_relative_error=match,
                fermion_count_rounded_unchanged=int(np.count_nonzero(ff==ff2)),
                exact_count_checks=nulls)


def run():
    P=[[.92,.22,.81],[.08,.78,.19]]
    F=[[.25,.71,.88],[.75,.29,.12]]
    W=[[.6,.2,.4],[.25,.25,.4],[.15,.55,.2]]
    a=[]
    for z in ([.57,-.63,.23,-.36],[-.7,.38,.31,-.53]):
        tests=[trial(e,.43,z,[.6,1.8,2.4],P,F,W) for e in (1e-3,1e-7,1e-12,1e-24)]
        v=tests[-1]['scaled_initial_signed_rates']
        if not min(v)<0<max(v) or tests[-1]['fermion_count_rounded_unchanged']!=4:
            raise ValueError(f"Opposing flux or sub-ulp chart not exercised: {z}, {v}.")
        a.append(dict(fermion_shifts=z,cases=tests))
    return dict(scope="Three manufactured shared Bose/Pauli events with normalized weights; no physical rates.",
                starts=a)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
