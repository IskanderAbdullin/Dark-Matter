"""An exact opt-in fast path for the endpoint capacity-limited collision step."""
import numpy as np

from early_pole_bounded_face_extents_pilot import step as reference


def step(ch,d,f,dt,*,rate_state=None):
    _,_,_,i,j,k,a,b,W=ch
    fs,fd,fp=np.split(f if rate_state is None else rate_state,(d.ns,d.iP))
    p,pc=d.parent_products(fp[k-1],fp[k],1.-b,b,False)
    q,qc=d.parent_products(fd[j-1],fd[j],1.-a,a,False)
    fo=W*p*(1.+fs[i])*qc
    inv=W*pc*fs[i]*q
    extent=dt*(fo-inv)
    ix=(i,d.ns+j-1,d.ns+j,d.iP+k-1,d.iP+k)
    dep=(np.ones_like(a),1-a,a,-(1-b),-b)
    flows=[x*extent for x in dep]
    out=sum(np.bincount(q,np.maximum(-s,0.),minlength=d.n) for q,s in zip(ix,flows,strict=True))
    inc=sum(np.bincount(q,np.maximum(s,0.),minlength=d.n) for q,s in zip(ix,flows,strict=True))
    stored=f*d.w
    holes=np.r_[np.full(d.ns,np.inf),d.w[d.ns:]*(1-f[d.ns:])]
    slack=1-32*np.finfo(float).eps
    # A sufficient componentwise condition for every original event
    # multiplier to clip to exactly one. At any other point retain the
    # original limiter, including its small numerical safety margin.
    if np.any(out>slack*stored) or np.any(inc>slack*holes):
        return reference(ch,d,f,dt,rate_state=rate_state)
    delta=sum(np.bincount(q,s,minlength=d.n) for q,s in zip(ix,flows,strict=True))
    new=f+delta/d.w
    if np.min(new)<-1e-16 or np.max(new[d.ns:])>1+1e-14:
        raise ArithmeticError('A capacity-proved step left the occupation bounds')
    return new,0.
