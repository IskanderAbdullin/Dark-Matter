"""True-event signed moment injection and convex realizability rejection."""
import json
import numpy as np
from scipy.special import expit
from true_event_fv_evolution_reference import energy_mesh, Reconstruction
from cell_number_energy_realizability_reference import fd_bounds


def check(M,n,mu_parent):
    q,m,b,(i,j,e,r,w)=energy_mesh(n,M,(.08,.18),6)
    B,F=Reconstruction(q,m,True),Reconstruction(q,b,False)
    gb=4/(2*np.pi**2)
    kb=e-m
    kf=r-b
    mu_b=-.02
    mu_f=.025
    sb,sf=mu_b-m,mu_f-b
    bf=1/np.expm1(B.k-sb)
    ff=expit(sf-F.k)
    nB=B.number*(B.w*bf).sum(1)
    nF=F.number*(F.w*ff).sum(1)
    kB=B.number*(B.w*bf*B.k).sum(1)
    kF=F.number*(F.w*ff*F.k).sum(1)
    p=expit(mu_parent-e-r)
    fB=1/np.expm1(e-mu_b)
    fF=expit(mu_f-r)
    J=w*(p*(1-fF)*(1+fB)-(1-p)*fF*fB)
    dNB=np.bincount(i,weights=J,minlength=n)
    dNF=np.bincount(j,weights=J,minlength=n)
    dKB=np.bincount(i,weights=kb*J,minlength=n)
    dKF=np.bincount(j,weights=kf*J,minlength=n)
    debit=-(e+r)@J
    assert abs(dNB.sum()-dNF.sum())<1e-14
    assert abs(dKB.sum()+dKF.sum()+m*dNB.sum()+b*dNF.sum()+debit)<1e-14
    def admissible(h):
        A, Bn=nB+h*dNB,nF+h*dNF
        Ka, Kb=kB+h*dKB,kF+h*dKF
        if np.any(A<0) or np.any(Bn<0) or np.any(Bn>F.number):
            return False,"number"
        if np.any(Ka<B.lo*A) or np.any(Ka>B.hi*A):
            return False,"boson_energy"
        for k,frac in enumerate(Bn/F.number):
            lo,hi=fd_bounds(q[k],q[k+1],b,frac)
            if not lo*gb<=Kb[k]<=hi*gb:
                return False,"fermion_energy"
        return True,"ok"
    assert admissible(0)[0]
    t,lo=1e-4,0.
    while admissible(t)[0]:
        lo=t
        t*=2
        assert t<1e12
    hi=t
    assert admissible(lo)[0]
    _,first=admissible(hi)
    for _ in range(40):
        mid=(lo+hi)/2
        if admissible(mid)[0]:
            lo=mid
        else:
            hi=mid
    assert admissible(lo)[0] and not admissible(hi)[0]
    limits=[]
    for nums,delta,upper in ((nB,dNB,None),(nF,dNF,F.number)):
        neg=delta<0
        if np.any(neg):
            limits.extend(nums[neg]/(-delta[neg]))
        pos=delta>0
        if upper is not None and np.any(pos):
            limits.extend((upper[pos]-nums[pos])/delta[pos])
    nlim=min(limits)
    assert lo<=nlim*(1+2e-12)
    return dict(M=M,n=n,mu_parent=mu_parent,number_only_step_limit=float(nlim),
                fraction_of_number_only_limit=float(lo/nlim),first_bound_reject=first,first_rejected_probe=t,
                largest_feasible_euler_trial=lo,smallest_rejected_euler_trial=hi,
                limiting_moment=admissible(hi)[1],
                relative_bracket_width=(hi-lo)/lo,
                absolute_number_source_error=float(abs(dNB.sum()-dNF.sum())),
                absolute_energy_source_error=float(abs(dKB.sum()+dKF.sum()+m*dNB.sum()+b*dNF.sum()+debit)))

if __name__ == "__main__":
    print(json.dumps({"cases":[check(m,n,u) for m in (.1,.005) for n in (12,24) for u in (.15,-.15)],
                      "status":"one signed state step; not a multistep evolution"},
                     indent=2,allow_nan=False))
