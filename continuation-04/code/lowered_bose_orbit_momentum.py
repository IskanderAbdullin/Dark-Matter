"""Momentum-constrained Onsager projection of the orbit/local event operator.
Entropy-positive at exact log mean, not eventwise kinematic conservation.
"""
import argparse,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from lowered_bose_orbit_collision import build

class Projector:
 def __init__(self,C,ng=8):
    self.C=C;self.ng=ng;self.p=C.orbit_p;self.shell=C.shell;self.ns=int(self.shell.max())+1
    assert np.all(self.shell[C.ids]==self.shell[C.ids[0]])
    self.localgroups=[np.flatnonzero(self.shell==s) for s in range(self.ns)]
    self.x,self.w=leggauss(ng);self.x=(self.x+1)/2;self.w/=2
 def project(self,f,diagnostic=False):
    C=self.C;M=C.A.shape[0];p=self.p;F,G,L=C.coeff(f)
    graw=((1+f)*(C.A.T@G)-f*(C.A.T@L))
    Kp=np.zeros((self.ns,len(f)));rawp=np.zeros(self.ns);raw=np.zeros_like(Kp)
    for s,ls in enumerate(self.localgroups):
        A=C.A[ls];raw[s]=(1+f)*(A.T@G[ls])-f*(A.T@L[ls]);rawp[s]=p@raw[s]
    lf=np.log(np.maximum(f,1e-300));lg=np.log1p(f)
    for u,wt in zip(self.x,self.w):
        xin=np.exp(u*lf+(1-u)*lg);xout=np.exp((1-u)*lf+u*lg)
        xs=[xin,xout];Z=[C.A@x for x in xs];P=[C.A@(x*p) for x in xs];mom=[np.divide(a,b,out=np.zeros_like(a),where=b>0) for a,b in zip(P,Z)]
        ids=C.ids;sg=np.array([-1.,-1.,1.,1.]);typ=[0,0,1,1];zz=np.array([Z[t][i] for t,i in zip(typ,ids)]);mm=np.array([mom[t][i] for t,i in zip(typ,ids)]);B=np.zeros((2,M));D=np.zeros((2,M))
        for j,t in enumerate(typ):
            others=[k for k in range(4) if k!=j];v=wt*C.weights*np.prod(zz[others],axis=0);h=sg[j]*np.sum(sg[others,None]*mm[others],axis=0)
            B[t]+=np.bincount(ids[j],weights=v,minlength=M);D[t]+=np.bincount(ids[j],weights=v*h,minlength=M)
        for s,ls in enumerate(self.localgroups):
            A=C.A[ls];Kp[s]+=sum(x*(p*(A.T@B[t,ls])+A.T@D[t,ls]) for t,x in enumerate(xs))
    den=Kp@p;lam=np.divide(rawp,den,out=np.zeros_like(den),where=den>1e-30);corr=lam@Kp;g=graw-corr
    if not diagnostic:return g/C.V
    ent=np.log1p(1/f);return g/C.V,dict(raw_pdot=float(p@graw),corrected_pdot=float(p@g),max_shell_pdot=float(abs(rawp-lam*den).max()),raw_entropy=float(ent@graw),corrected_entropy=float(ent@g),entropy_removed=float(ent@corr),correction_L1=float(abs(corr).sum()),raw_L1=float(abs(graw).sum()),max_lambda=float(abs(lam).max()),min_den=float(den.min()),number=float(g.sum()),energy=float(C.e@g),Kp_number=float(abs(Kp.sum(1)).max()),Kp_energy=float(abs(Kp@C.e).max()),raw_decomposition=float(abs(raw.sum(0)-graw).max()),entropy_formula_error=float(abs(ent@corr - np.sum(lam*rawp))))

def run(NE,K,nL,nphase,m,radsub,ng):
    C,diag=build(NE,K,nL,nphase,m,radsub);P=Projector(C,ng);be=1/np.expm1(.5+C.e-.1*C.orbit_p);g,b=P.project(be,True);raw=C.derivative(be);diag.update(moving_BE=dict(raw_max=float(abs(raw).max()),corrected_max=float(abs(g).max()),**b))
    f=1/np.expm1(.5+C.e)*(1+.3*np.sin(2*np.pi*C.phase))*(1+.15*np.cos(2*np.pi*C.phase));g,b=P.project(f,True);diag['phase_perturbation']=b
    rng=np.random.default_rng(51);f=1/np.expm1(.5+C.e)*np.exp(rng.uniform(-20,0,len(C.e)));g,b=P.project(f,True);diag['low_occupation']=dict(min=float(f.min()),neg_rate_below_1e_10=int(np.sum((f<1e-10)&(g<0))),min_forward_euler_dt=float(np.min(-f[g<0]/g[g<0])),**b)
    return diag
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=6);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=8);p.add_argument('--nphase',type=int,default=32);p.add_argument('--m',type=int,default=7);p.add_argument('--radsub',type=int,default=1);p.add_argument('--ng',type=int,default=8);a=p.parse_args();print(json.dumps(run(a.NE,a.K,a.nL,a.nphase,a.m,a.radsub,a.ng)))
