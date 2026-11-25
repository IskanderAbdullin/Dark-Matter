"""Bound harmonic-trap demonstration of orbit/local conservative nonlinear UU coupling.
Not a self-gravitating halo or unbound solver. All units dimensionless.
"""
import argparse,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix,diags
from scipy.stats import qmc
from lowered_bose_projected_axis import hats,rounded_distribution

def shift(f,a,order=2):
    k=int(np.floor(a));x=a-k;y=np.roll(f,k);prev=np.roll(y,1);l=y-prev;r=np.roll(y,-1)-y
    slope=np.where(l*r>0,np.sign(l)*np.minimum(abs(l),abs(r)),0.) if order==2 else np.zeros_like(y)
    return (1-x)*y+x*prev+x*(1-x)/2*(np.roll(slope,1)-slope)

class Coupling:
 def __init__(self,V,e,W,local_e,p,ids,weights,orbit_p=None):
    self.V=np.asarray(V);self.e=e;self.W=W.tocsr();self.local_e=local_e;self.p=p;self.ids=ids;self.weights=weights
    self.local_mass=np.asarray(W.sum(1)).ravel();self.orbit_p=orbit_p
    self.A=diags(np.divide(1,self.local_mass,out=np.zeros_like(self.local_mass),where=self.local_mass>0))@self.W
    assert np.allclose(np.asarray(W.sum(0)).ravel(),V,rtol=1e-12,atol=1e-14)
    assert np.max(np.abs(e[W.nonzero()[1]]-local_e[W.nonzero()[0]]))<1e-12
    assert np.max(abs(local_e[ids[0]]+local_e[ids[1]]-local_e[ids[2]]-local_e[ids[3]]))<1e-12
    assert np.max(abs(p[ids[0]]+p[ids[1]]-p[ids[2]]-p[ids[3]]))<1e-12
 def coeff(self,f):
    F=self.A@f;a,b,c,d=F[self.ids];w=self.weights;M=len(F);G=np.zeros(M);L=np.zeros(M)
    for j,(g,l) in enumerate([(c*d*(1+b),b*(1+c)*(1+d)),(c*d*(1+a),a*(1+c)*(1+d)),(a*b*(1+d),d*(1+a)*(1+b)),(a*b*(1+c),c*(1+a)*(1+b))]):
        G+=np.bincount(self.ids[j],weights=w*g,minlength=M);L+=np.bincount(self.ids[j],weights=w*l,minlength=M)
    return F,G,L
 def derivative(self,f,diagnostic=False):
    F,G,L=self.coeff(f);gain=self.A.T@G;loss=self.A.T@L
    g=((1+f)*gain-f*loss)/self.V
    if not diagnostic:return g
    local=(1+F)*G-F*L
    transfer=self.A.multiply((1+f)[None,:]).multiply(G[:,None])-self.A.multiply(f[None,:]).multiply(L[:,None])
    assert np.allclose(np.asarray(transfer.sum(1)).ravel(),local,rtol=1e-10,atol=1e-12)
    assert np.allclose(np.asarray(transfer.sum(0)).ravel(),g*self.V,rtol=1e-10,atol=1e-12)
    return g,dict(local_ledger=[float(local.sum()),float(self.local_e@local),float(self.p@local)],physical_pdot=float((self.V*self.orbit_p)@g) if self.orbit_p is not None else None,absolute_local_source=float(abs(local).sum()))
 def entropy(self,f):return float(self.V@((1+f)*np.log1p(f)-f*np.log(np.maximum(f,1e-300))))
 def stats(self,f):return dict(N=float(self.V@f),E=float((self.V*self.e)@f),S=self.entropy(f),min=float(f.min()),p=float((self.V*self.orbit_p)@f) if self.orbit_p is not None else None)

def brute_control(seed=42):
    rng=np.random.default_rng(seed);W=np.zeros((4,8));W[np.arange(4)[:,None],np.arange(8).reshape(4,2)]=rng.uniform(.2,1,(4,2));e=np.repeat([1.,4.,2.,3.],2);p=np.array([-1.,1.,-.5,.5]);V=W.sum(0);ids=np.arange(4)[:,None];C=Coupling(V,e,coo_matrix(W),e[::2],p,ids,np.array([.7]));f=rng.lognormal(0,1,8);g=C.derivative(f);direct=np.zeros(8);S=0.
    for i in np.ndindex(2,2,2,2):
        ix=2*np.arange(4)+i;wt=.7*np.prod(W[np.arange(4),ix]/W.sum(1));a,b,c,d=f[ix];X=a*b*(1+c)*(1+d);Y=c*d*(1+a)*(1+b);j=wt*(X-Y);direct[ix]+=j*np.array([-1,-1,1,1]);S+=j*np.log(X/Y)
    return dict(max_source_error=float(abs(direct-g*V).max()),entropy_rate=float((V*np.log1p(1/f))@g),entropy_event=S,number=float(V@g),energy=float((V*e)@g))

def build(NE,K,nL,nphase,m=7,radsub=1,Emax=6.,seed=0):
    de=Emax/NE;en=(np.arange(NE)+.5)*de;dp=np.sqrt(Emax)/K;xe,we=leggauss(nL);x,w=leggauss(6);T=np.pi
    edges=np.r_[0,en];U=np.concatenate([np.linspace(a,b,radsub+1)[:-1] for a,b in zip(edges[:-1],edges[1:])]+[edges[-1:]])
    R=np.sqrt(2*U);NR=len(R)-1;tab=np.full((NR,NE,2*K+1),-1,int);le=[];lp=[];expected=[];rc=[];local_ranges=[]
    for j,(lo,hi) in enumerate(zip(U[:-1],U[1:])):
        u=(lo+hi)/2;vr=(R[j+1]**3-R[j]**3)/3;active=np.flatnonzero(en>=hi-1e-12);start=len(le)
        for i in active:
            E=en[i]-u;mom=hats(np.sqrt(E),dp,K);sel=np.flatnonzero(mom[:,0]>1e-16);tab[j,i,sel]=np.arange(len(le),len(le)+len(sel));le.extend([en[i]]*len(sel));lp.extend((sel-K)*dp);expected.extend(vr*de*mom[sel,0]/4);rc.extend([j]*len(sel))
        local_ranges.append((active,u,vr,start,len(le)))
    le=np.array(le);lp=np.array(lp);expected=np.array(expected);rc=np.array(rc);nlocal=len(le);V=[];oe=[];op=[];oL=[];orph=[];rows=[];cols=[];data=[];norb=NE*nL;phase=np.linspace(0,1,nphase+1);mid=(phase[:-1]+phase[1:])/2
    for i,e in enumerate(en):
        for z,ww in zip(xe,we):
            L=(z+1)*e*e/2;WL=ww*e*e/2;A=np.sqrt(e*e-L);base=len(V);volume=de*WL*T/(8*np.sqrt(2));V.extend([volume/nphase]*nphase);oe.extend([e]*nphase);oL.extend([L]*nphase);orph.extend(mid)
            th=mid[:,None]+x[None,:]/(2*nphase);r=np.sqrt(e-A*np.cos(2*np.pi*th));p=A*np.sin(2*np.pi*th)/(np.sqrt(2)*r);op.extend(p@w/2)
            for j,(lo,hi) in enumerate(zip(U[:-1],U[1:])):
                if hi>e+1e-12:continue
                l=np.arccos(np.clip((e-2*lo)/A,-1,1))/(2*np.pi);h=np.arccos(np.clip((e-2*hi)/A,-1,1))/(2*np.pi)
                if h<=l:continue
                u=(lo+hi)/2;E=en[i]-u
                for a,b in [(l,h),(1-h,1-l)]:
                    ii=np.arange(max(0,int(a*nphase)),min(nphase,int(np.ceil(b*nphase))));aa=np.maximum(a,phase[ii]);bb=np.minimum(b,phase[ii+1]);keep=bb>aa;ii=ii[keep];aa=aa[keep];bb=bb[keep]
                    th=(aa[:,None]+bb[:,None])/2+(bb-aa)[:,None]*x/2;r=np.sqrt(e-A*np.cos(2*np.pi*th));p=A*np.sin(2*np.pi*th)/(np.sqrt(2)*r);q=np.sqrt(np.maximum(e-r*r/2,1e-300));val=np.clip(p/q,-1,1)*np.sqrt(E)/dp;ix=np.floor(val).astype(int);f=val-ix;wt=volume*(bb-aa)[:,None]*w/2
                    for k,zwt in [(ix,1-f),(ix+1,f)]:
                        mask=zwt>1e-15;loc=tab[j,i,k[mask]+K];assert np.all(loc>=0);rows.extend(loc);cols.extend(np.broadcast_to((base+ii)[:,None],th.shape)[mask]);data.extend((wt*zwt)[mask])
    V=np.array(V);oe=np.array(oe);op=np.array(op);oL=np.array(oL);orph=np.array(orph);W=coo_matrix((data,(rows,cols)),shape=(nlocal,len(V))).tocsr();local_mass=np.asarray(W.sum(1)).ravel();missing=local_mass<=1e-20
    samples=qmc.Sobol(4,scramble=True,seed=seed).random_base2(m);allids=[];allw=[];lost=0.;total=0.
    for j,(active,u,vr,_,_) in enumerate(local_ranges):
        n=len(active)
        for s in range(2*n-1):
            pairs=np.arange(max(0,s-n+1),min(n-1,s//2)+1);other=s-pairs
            for a,b in zip(*np.triu_indices(len(pairs),0)):
                es=active[np.array([pairs[a],other[a],pairs[b],other[b]])];E=en[es]-u;momidx,v=rounded_distribution(E,dp,K,samples);ids=tab[j,es[:,None],momidx+K];assert np.all(ids>=0)
                wq=vr*de**3*np.sqrt(E).min()/8*(1 if es[0]==es[1] else 2)*(1 if es[2]==es[3] else 2)*(.5 if a==b else 1.)*v/v.sum();keep=~np.any(missing[ids],axis=0);lost+=wq[~keep].sum();total+=wq.sum();allids.append(ids[:,keep]);allw.append(wq[keep])
    C=Coupling(V,oe,W,le,lp,np.concatenate(allids,axis=1),np.concatenate(allw),op)
    C.nphase=nphase;C.norb=norb;C.oL=oL;C.phase=orph
    pp=(W.T@lp)/V
    diag=dict(NE=NE,K=K,nL=nL,nphase=nphase,m=m,radsub=radsub,states=len(V),local_states=nlocal,events=len(C.weights),missing_local=int(missing.sum()),lost_event_weight=float(lost/total),local_mass_l1=float(abs(local_mass-expected).sum()/expected.sum()),total_volume_rel=float(V.sum()/expected.sum()-1),p_map_rms=float(np.sqrt((V@(pp-op)**2)/V.sum())),overlap_column_relative=float(np.max(abs(np.asarray(W.sum(0)).ravel()/V-1))))
    return C,diag

def run(NE,K,nL,nphase,m,radsub,tmax,dt):
    C,diag=build(NE,K,nL,nphase,m,radsub);f0=1/np.expm1(.5+C.e);eq=C.derivative(f0);f=f0*(1+.3*np.sin(2*np.pi*C.phase))*(1+.15*np.cos(2*np.pi*C.phase));g,a=C.derivative(f,True);diag.update(eq_max=float(abs(eq).max()),initial=C.stats(f),initial_collision=a,entropy_rate=float((C.V*np.log1p(1/f))@g),brute=brute_control())
    t=0;steps=0;reject=0;min_dS=0;max_mom=0
    def stream(f,t):return np.concatenate([shift(row,t*nphase/np.pi,2) for row in f.reshape(-1,nphase)])
    fs=f.copy()
    while t<tmax-1e-12:
        d=min(dt,tmax-t)
        while True:
            x=stream(f,d/2);g=C.derivative(x);p=x+d*g
            if p.min()>=0:
                y=.5*(x+p+d*C.derivative(p));ds=C.entropy(y)-C.entropy(x) if y.min()>=0 else -1.
                if y.min()>=0 and ds>=-1e-12:break
            d*=.5;reject+=1
        min_dS=min(min_dS,ds);max_mom=max(max_mom,abs((C.V*C.orbit_p)@(y-x)));f=stream(y,d/2);fs=stream(fs,d);steps+=1;t+=d
    diag.update(final=C.stats(f),pure_stream=C.stats(fs),difference_l1=float(C.V@abs(f-fs)),min_collision_entropy_diff=min_dS,max_collision_physical_p_increment=float(max_mom),t=t,steps=steps,rejects=reject)
    return diag
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--NE',type=int,default=6);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=8);p.add_argument('--nphase',type=int,default=32);p.add_argument('--m',type=int,default=7);p.add_argument('--radsub',type=int,default=1);p.add_argument('--t',type=float,default=.2);p.add_argument('--dt',type=float,default=.02);a=p.parse_args();print(json.dumps(run(a.NE,a.K,a.nL,a.nphase,a.m,a.radsub,a.t,a.dt)))
