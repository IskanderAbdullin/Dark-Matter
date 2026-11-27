"""Pilot: resolved nonlinear one-event limit when fermion count changes are sub-ulp."""
import json

import numpy as np
from scipy.integrate import quad,solve_ivp

DT=np.longdouble


def trial(e,p,q,H=.75):
    eps, P, Q = map(DT,(e,p,q))
    f=DT(.43)
    u=f*(1-f)
    def sd(x):
        x=DT(x)
        parent=P-x; final=Q+x
        def perturb(d):
            return np.log1p(-eps*d/(1-f))-np.log1p(eps*d/f)
        # n_B(ξ)=1/eps+eps*x, with ξ=eps*x and time=eps*s.
        y=np.log1p(eps/(1+eps**2*x))
        A=y+perturb(final)-perturb(parent)
        ph=(f+eps*parent)*(1-f-eps*final)
        if ph<=0:
            raise ValueError("Toy perturbation crossed a Fermi boundary.")
        ahead=ph*eps/(-np.expm1(-y))
        S=ahead*(-np.expm1(-A))/eps
        return float(S),float(A/eps)

    v=solve_ivp(lambda t,z:[sd(z[0])[0]],(0,H),[0.],method="DOP853",rtol=9e-13,atol=3e-15)
    if not v.success: raise RuntimeError(v.message)
    x=v.y[0,-1]
    expected=(u+P-Q)*(-np.expm1(-2*DT(H)))/2
    es=quad(lambda s:sd(s)[1],0,x,epsabs=1e-12,epsrel=1e-12)[0]
    if es<=0 or np.sign(x)!=np.sign(expected):
        raise ValueError("Perturbative extent or scaled entropy sign failed.")
    import mpmath as mp
    def C(x):
        n,d=DT(x).as_integer_ratio()
        return mp.mpf(n)/d
    with mp.workdps(120):
        ee,vv,pp,qq,xx=map(C,(eps,f,P,Q,x))
        L=lambda n:-n*mp.log(n)-(1-n)*mp.log1p(-n)
        Qs=lambda n:mp.log1p(n)+n*mp.log1p(1/n)
        before=[vv+ee*pp,vv+ee*qq,1/ee]
        after=[before[0]-ee*xx,before[1]+ee*xx,before[2]+ee*xx]
        check=(sum(L(y)-L(z) for y,z in zip(after[:2],before[:2],strict=True))+Qs(after[2])-Qs(before[2]))/ee**2
    err=abs(es/float(check)-1)
    if err>2e-10: raise ValueError("Independent high-precision scaled entropy did not close.")
    if e<=1e-24 and (abs(x-expected)>3e-13 or float(f+eps*(P-x))!=float(f+eps*P)):
        raise ValueError("Sub-ulp count fixture or independent leading limit failed.")
    return dict(inverse_bose_count=e,scaled_extent=float(x),leading_solution=float(expected),
                leading_error=float(abs(x-expected)),entropy_increase_over_inverse_count_squared=es,
                independent_scaled_entropy_error=err,
                naive_rounded_fermi_change=float(f+eps*(P-x))==float(f+eps*P))


if __name__=='__main__':
    print(json.dumps(dict(scope="Manufactured nonlinear single-event scaled limit; no physical abundance.",
         cases=[dict(start=[p,q],rows=[trial(e,p,q) for e in (1e-3,1e-7,1e-12,1e-24)])
                for p,q in ((.34,-.1),(-.5,.1))]),indent=2,allow_nan=False))
