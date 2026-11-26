"""Same projected prompt model with parent-local assembled Jacobian."""
import numpy as np
from joint_resonant_centered_stiff_reference import Centered


def compile_prompt(model):
    k=model.k
    segments=[]
    if not len(k.y): return segments
    ends=np.r_[0,np.flatnonzero(np.diff(k.y))+1,len(k.y)]
    for l,r in zip(ends,ends[1:]):
        s=slice(l,r)
        idx=np.stack([2*k.I[s],2*k.I[s]+1,2*(model.n+k.J[s]),2*(model.n+k.J[s])+1],axis=1)
        ix,inv=np.unique(idx,return_inverse=True)
        flow=np.empty(idx.shape)
        flow[:,0]=1.;flow[:,1]=model.us[s]
        flow[:,2]=1.;flow[:,3]=model.uc[s]
        P=np.zeros((r-l,len(ix)))
        P[np.arange(r-l)[:,None],inv.reshape(idx.shape)]=flow
        segments.append((s,ix,P,ix<2*model.n))
    return segments


def local_prompt(model,z,qref,R0,segments):
    k=model.k
    f,g=model.arrays(qref+z)
    f=k.c[0].f(f,k.ph,k.I);g=k.c[1].f(g,k.ch,k.J)
    b=k.pre*(1+f)*(1-g)
    x,y=model.arrays(z)
    shift=x[0,k.I]+x[1,k.I]*model.us+y[0,k.J]+y[1,k.J]*model.uc
    md=np.expm1(-shift)
    # Keep the value path identical to the established operator.
    def accum(v):return np.bincount(k.y,weights=v,minlength=k.n)[k.y]
    B=accum(b);J=accum(b*md);D=B+R0*(B+J)
    omega=np.divide(k.fac*R0,D*model.H,out=np.zeros_like(D),where=D>0.)
    v=J-B*md
    out=model.proj(omega*b*v)
    cache=None
    def derivative(d):
        nonlocal cache
        if cache is None:
            aa=np.zeros((model.d,model.d))
            for s,ix,P,sf in segments:
                m=md[s,None];w=b[s,None]; vv=v[s,None]
                db=w*(P*np.where(sf,-f[s,None],g[s,None]))
                dm=-(1+m)*P
                dB=db.sum(axis=0)
                dJ=(db*m+w*dm).sum(axis=0)
                dD=dB+R0[s.start]*(dB+dJ)
                change=dJ-m*dB-B[s,None]*dm
                sink=np.divide(dD,D[s.start],out=np.zeros_like(dD),where=D[s.start]>0.)
                df=omega[s,None]*(vv*db+w*change-w*vv*sink)
                aa[np.ix_(ix,ix)]+=P.T@df
            cache=aa
        return cache@d
    return out,derivative


class LocalPrompt(Centered):
    def centered(self,z,qref,R0):
        if not hasattr(self,"prompt_segments"):
            self.prompt_segments=compile_prompt(self)
        return local_prompt(self,z,qref,R0,self.prompt_segments)
