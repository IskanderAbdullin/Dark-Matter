"""Batched conservative frozen-orbit phase remap. Geometry cached by dt."""
import numpy as np

class Batch:
 def __init__(self,edges,open_mask,T,V,e,phase_mass):
    self.edges=np.asarray(edges);self.open=bool(open_mask);self.T=np.asarray(T);self.V=np.asarray(V);self.e=np.asarray(e);self.phase_mass=np.asarray(phase_mass);self.d=np.diff(edges);self.x=(edges[1:]+edges[:-1])/2;self.left=np.r_[self.x[-1]-1,self.x[:-1]];self.right=np.r_[self.x[1:],self.x[0]+1];self.cache={}
 def geometry(self,dt):
    if dt in self.cache:return self.cache[dt]
    a=dt/self.T[:,None];a=a if self.open else a%1
    lo=self.edges[:-1][None,:]-a;hi=self.edges[1:][None,:]-a
    if self.open:
        lo=np.clip(lo,0,1);hi=np.clip(hi,0,1);nturn=np.zeros_like(lo);same_turn=np.ones_like(lo,dtype=bool)
    else:
        tl=np.floor(lo);th=np.floor(hi);lo=lo-tl;hi=hi-th;nturn=th-tl;same_turn=tl==th
    i=np.clip(np.searchsorted(self.edges,lo,side='right')-1,0,len(self.d)-1);j=np.clip(np.searchsorted(self.edges,hi,side='right')-1,0,len(self.d)-1);A=lo-self.edges[i];B=hi-self.edges[j];same=(i==j)&same_turn
    val=(i,j,A,B,nturn,same)
    if len(self.cache)>8:self.cache.clear()
    self.cache[dt]=val;return val
 def apply(self,f,dt,order=2):
    d=self.d;x=self.x;g=np.roll(f,1,axis=1);h=np.roll(f,-1,axis=1);a=(f-g)/(x-self.left);b=(h-f)/(self.right-x);s=np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.) if order==2 else np.zeros_like(f);s=np.clip(s,-2*f/d,2*f/d)
    if self.open:s[:,0]=s[:,-1]=0.
    pre=np.concatenate((np.zeros((len(f),1)),np.cumsum(f*d,axis=1)),axis=1);i,j,A,B,nturn,same=self.geometry(dt);get=lambda x,ind:np.take_along_axis(x,ind,axis=1);fi=get(f,i);fj=get(f,j);si=get(s,i);sj=get(s,j);di=d[i];dj=d[j]
    v=(di-A)*(fi+si*A/2)+B*(fj+sj*(B-dj)/2)+get(pre,j)-get(pre,i+1)+nturn*pre[:,-1,None]
    v=np.where(same,(B-A)*(fi+si*(A+B-di)/2),v)
    y=v/d;out=self.phase_mass*(pre[:,-1]-v.sum(1)) if self.open else np.zeros(len(f));return y,out

class Streamer:
 def __init__(self,C,TIME):
    self.C=C;self.TIME=TIME;self.n=C.nphase;self.groups=[]
    for isopen in [False,True]:
        ii=np.array([i for i,(o,_,_) in enumerate(C.orbits) if o.open==isopen],int)
        if not len(ii):continue
        edge=C.orbits[ii[0]][1];assert all(np.array_equal(C.orbits[i][1],edge) for i in ii);T=np.array([C.orbits[i][0].T for i in ii]);ee=np.array([-C.orbits[i][0].s for i in ii]);vv=C.omass[ii]
        self.groups.append((ii,Batch(edge,isopen,T,None,ee,vv)))
 def apply(self,f,dt,order=2):
    n=self.n;z=f.reshape(-1,n);g=np.empty_like(z);out=np.zeros(3);mn=0.
    for ii,B in self.groups:
        y,w=B.apply(z[ii],dt/self.TIME,order);g[ii]=y;out+=np.array([w.sum(),w@B.e,w@(B.e>self.C.e0)]);mn=min(mn,float(y.min()))
    if mn<-1e-9:raise RuntimeError('negative batched stream')
    return np.maximum(g.ravel(),0),out,mn
