"""Positive reversible-flux tilt enforcing physical orbit-cell p per shell.
Not eventwise momentum conservation; a constrained numerical closure.
"""
import argparse,json
import numpy as np
from scipy.optimize import brentq
from lowered_bose_orbit_collision import build,shift

class Tilt:
 def __init__(self,C):
    self.C=C;self.p=C.orbit_p;self.groups=[]
    for s in range(int(C.shell.max())+1):
        ls=np.flatnonzero(C.shell==s);es=np.flatnonzero(C.shell[C.ids[0]]==s);assert np.all(C.shell[C.ids[:,es]]==s)
        if not len(es):continue
        A=C.A[ls].tocsr();ix=np.unique(A.indices);A=A[:,ix].tocsr();ids=np.searchsorted(ls,C.ids[:,es]);self.groups.append((A,ix,ids,C.weights[es]))
 def evaluate(self,A,ix,ids,w,f,lam,source=False):
    p=self.p[ix];a=f[ix]*np.exp(lam*p/2);b=(1+f[ix])*np.exp(-lam*p/2);F=A@a;B=A@b
    if not (np.isfinite(F).all() and np.isfinite(B).all()):raise FloatingPointError('tilt overflow')
    u=np.array([F[ids[0]],F[ids[1]],B[ids[2]],B[ids[3]]]);v=np.array([B[ids[0]],B[ids[1]],F[ids[2]],F[ids[3]]]);X=np.prod(u,axis=0);Y=np.prod(v,axis=0)
    if source:
        gain=np.zeros(A.shape[0]);loss=np.zeros_like(gain)
        for j in range(4):
            k=[i for i in range(4) if i!=j];g=np.prod(v[k],axis=0) if j<2 else np.prod(u[k],axis=0);l=np.prod(u[k],axis=0) if j<2 else np.prod(v[k],axis=0)
            gain+=np.bincount(ids[j],weights=w*g,minlength=len(gain));loss+=np.bincount(ids[j],weights=w*l,minlength=len(loss))
        return b*(A.T@gain)-a*(A.T@loss)
    P=A@(a*p);Q=A@(b*p);m=np.divide(P,F,out=np.zeros_like(P),where=F>0);n=np.divide(Q,B,out=np.zeros_like(Q),where=B>0)
    dx=-m[ids[0]]-m[ids[1]]+n[ids[2]]+n[ids[3]];dy=-n[ids[0]]-n[ids[1]]+m[ids[2]]+m[ids[3]]
    return float(w@(X*dx-Y*dy))
 def derivative(self,f,diagnostic=False):
    C=self.C;g=np.zeros_like(f);raw=np.zeros_like(f);maxp=0.;lams=[];rawp=[];evals=0
    for A,ix,ids,w in self.groups:
        fn=lambda l:self.evaluate(A,ix,ids,w,f,l)
        v=fn(0.);rawp.append(v);evals+=1
        if abs(v)<1e-18:lam=0.
        else:
            b=1. if v>0 else -1.;vb=fn(b);evals+=1
            while np.sign(vb)==np.sign(v) and abs(b)<256:
                b*=2;vb=fn(b);evals+=1
            if np.sign(vb)==np.sign(v):raise RuntimeError('no finite momentum root on bracket; zero-support or stiff constraint')
            lam,res=brentq(fn,min(0,b),max(0,b),xtol=1e-12,full_output=True);evals+=res.function_calls
        z=self.evaluate(A,ix,ids,w,f,lam,True);g[ix]+=z;maxp=max(maxp,abs(self.p[ix]@z));lams.append(lam)
        if diagnostic:raw[ix]+=self.evaluate(A,ix,ids,w,f,0.,True)
    if not diagnostic:return g/C.V
    ent=np.log1p(1/np.maximum(f,1e-300));return g/C.V,dict(max_lambda=float(max(abs(np.array(lams)))),max_shell_pdot=float(maxp),raw_pdot=float(self.p@raw),corrected_pdot=float(self.p@g),raw_entropy=float(ent@raw),corrected_entropy=float(ent@g),correction_L1=float(abs(g-raw).sum()),raw_L1=float(abs(raw).sum()),number=float(g.sum()),energy=float(C.e@g),raw_reference=float(abs(raw-C.derivative(f)*C.V).max()),root_evals=evals)

def run(NE,K,nL,nphase,m,radsub,tmax,dt):
    C,diag=build(NE,K,nL,nphase,m,radsub);P=Tilt(C);be=1/np.expm1(.5+C.e-.1*C.orbit_p);g,b=P.derivative(be,True);diag.update(moving_BE=dict(raw_max=float(abs(C.derivative(be)).max()),corrected_max=float(abs(g).max()),**b))
    f=1/np.expm1(.5+C.e)*(1+.3*np.sin(2*np.pi*C.phase))*(1+.15*np.cos(2*np.pi*C.phase));g,b=P.derivative(f,True);diag['phase_perturbation']=b
    rng=np.random.default_rng(51);z=1/np.expm1(.5+C.e)*np.exp(rng.uniform(-20,0,len(C.e)));z[rng.random(len(z))<.2]=0.;g,b=P.derivative(z,True);diag['zero_occupation']=dict(zeros=int(np.sum(z==0)),min_zero_rate=float(g[z==0].min()),min_forward_euler_dt=float(np.min(-z[g<0]/g[g<0])),**b)
    initial=C.stats(f);t=0.;steps=0;reject=0;min_ds=0;max_p=0.;fs=f.copy();fr=f.copy()
    def stream(f,t):return np.concatenate([shift(row,t*nphase/np.pi,2) for row in f.reshape(-1,nphase)])
    while t<tmax-1e-12:
        d=min(dt,tmax-t)
        while True:
            x=stream(f,d/2);g=P.derivative(x);p=x+d*g
            if p.min()>=0:
                y=.5*(x+p+d*P.derivative(p));ds=C.entropy(y)-C.entropy(x) if y.min()>=0 else -1.
                if y.min()>=0 and ds>=-1e-12:break
            d*=.5;reject+=1
        min_ds=min(min_ds,ds);max_p=max(max_p,abs((C.V*C.orbit_p)@(y-x)));f=stream(y,d/2);fs=stream(stream(fs,d/2),d/2)
        x=stream(fr,d/2);p=x+d*C.derivative(x);fr=stream(.5*(x+p+d*C.derivative(p)),d/2)
        steps+=1;t+=d
    diag.update(initial=initial,final=C.stats(f),pure_stream=C.stats(fs),raw_final=C.stats(fr),physical_change_L1=float(C.V@abs(f-fs)),correction_effect_L1=float(C.V@abs(f-fr)),max_collision_p_increment=float(max_p),min_collision_entropy_diff=min_ds,t=t,steps=steps,rejects=reject)
    return diag
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=6);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=8);p.add_argument('--nphase',type=int,default=32);p.add_argument('--m',type=int,default=7);p.add_argument('--radsub',type=int,default=1);p.add_argument('--t',type=float,default=.2);p.add_argument('--dt',type=float,default=.02);a=p.parse_args();print(json.dumps(run(a.NE,a.K,a.nL,a.nphase,a.m,a.radsub,a.t,a.dt)))
