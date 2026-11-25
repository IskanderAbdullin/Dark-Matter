"""Frozen self-consistent *initial* lowered Bose potential, open/closed orbit transport.
Finite domain/energy cutoffs; no evolving Poisson. Conditional discretization study.
"""
import argparse,json,math,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq,minimize_scalar
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix
from scipy.stats import qmc
from lowered_antonov_hessian import Halo
from lowered_bose_projected_axis import hats,rounded_distribution
from lowered_bose_orbit_collision import Coupling
from lowered_bose_orbit_tilt import Tilt
from lowered_bose_orbit_remap import mesh,shift,integrator,slopes
R=.6*3.085677581491367e19;V=44100.;GYR=1e9*365.25*86400;TIME=R/V/GYR

def potential(H):
    e0=H.Qt/H.rt
    def psi(r):
        r=np.asarray(r);return np.where(r<=H.rt,H.psi(np.minimum(r,H.rt).ravel()).reshape(r.shape),H.Qt/np.maximum(r,H.rt)-e0)
    return psi

def occ(H,s):return H.F(np.maximum(s,0))

class Orbit:
 def __init__(self,psi,s,L,rmax,rc,isopen,nz=128):
    self.s=s;self.L=L;self.open=isopen;self.psi=psi
    g=lambda r:2*r*r*(float(psi(r))-s)-L
    lo=brentq(g,0,rc,xtol=1e-12);hi=rmax if isopen else brentq(g,rc,rmax,xtol=1e-12);self.lo=lo;self.hi=hi
    z=np.linspace(0,np.pi/2,nz+1);x,w=leggauss(6);a=(z[:-1,None]+z[1:,None])/2+np.diff(z)[:,None]*x/2;rr=lo+(hi-lo)*np.sin(a)**2
    vr=np.sqrt(np.maximum(2*(psi(rr)-s)-L/rr**2,1e-300));d=(hi-lo)*np.sin(2*a)/vr;tt=np.r_[0,np.cumsum(np.diff(z)/2*(d@w))];self.T=2*tt[-1]
    assert np.isfinite(self.T) and self.T>0 and np.all(np.diff(tt)>0)
    self.inv=PchipInterpolator(tt,lo+(hi-lo)*np.sin(z)**2);self.forward=PchipInterpolator(lo+(hi-lo)*np.sin(z)**2,tt)
 def coords(self,th):
    th=np.asarray(th);t=self.T*(abs(th-.5) if self.open else np.minimum(th,1-th));r=self.inv(t);p=np.sign(th-.5 if self.open else .5-th)*np.sqrt(np.maximum(self.psi(r)-self.s-self.L/(2*r*r),0));return r,p
 def intervals(self,a,b):
    a=max(a,self.lo);b=min(b,self.hi)
    if b<=a:return []
    l=float(self.forward(a)/self.T);h=float(self.forward(b)/self.T)
    return [(max(0,.5-h),min(.5,.5-l)),(max(.5,.5+l),min(1,.5+h))] if self.open else [(max(0,l),min(.5,h)),(max(.5,1-h),min(1,1-l))]

def open_shift(edges,f,dt,order=2):
    d=np.diff(edges);s=slopes(edges,f,order);s[0]=s[-1]=0.;prefix=np.r_[0,np.cumsum(d*f)]
    lo=np.clip(edges[:-1]-dt,0,1);hi=np.clip(edges[1:]-dt,0,1);i=np.clip(np.searchsorted(edges,lo,side='right')-1,0,len(f)-1);j=np.clip(np.searchsorted(edges,hi,side='right')-1,0,len(f)-1);a=lo-edges[i];b=hi-edges[j];same=i==j
    val=(d[i]-a)*(f[i]+s[i]*a/2)+b*(f[j]+s[j]*(b-d[j])/2)+prefix[j]-prefix[i+1];val[same]=(b[same]-a[same])*(f[i[same]]+s[i[same]]*(a[same]+b[same]-d[same])/2)
    g=val/d;out=float(prefix[-1]-val.sum());return g,out

def reference(H,n=24,nr=96):
    """Independent continuum initial above-cutoff and true-unbound gains."""
    x,w=leggauss(n);x=(x+1)/2;w/=2;u,ww=leggauss(nr);u=(u+1)/2;ww/=2;r=H.rt*u*u;vr=2*H.rt*u*r*r*ww;psi=H.psi(r);e0=H.Qt/H.rt;rows=[]
    for ex in [0.,e0]:
        total=0.
        for ps,radw in zip(psi,vr):
            if ps<=ex:continue
            gap=H.eta-ps;la=np.log((gap+ps)/(gap+ex));a=(gap+ex)*np.exp(la*x[:,None,None])-gap;wa=(gap+a)*la*w[:,None,None]
            blo=ps+ex-a;lb=np.log((gap+ps)/(gap+blo));b=(gap+blo)*np.exp(lb*x[None,:,None])-gap;wb=(gap+b)*lb*w[None,:,None]
            cmax=a+b-ps-ex;lc=np.log1p(cmax/gap);c=gap*np.expm1(lc*x[None,None,:]);wc=(gap+c)*lc*w[None,None,:]
            f1=occ(H,ps-a);f2=occ(H,ps-b);f3=occ(H,ps-c);total+=radw*np.sum(wa*wb*wc*np.sqrt(np.minimum(np.minimum(a,b),c))*f1*f2*(1+f3))/4
        rows.append(float(total))
    return np.array(rows)

def build(H,M=1,K=4,nL=6,nphase=64,nrad=12,m=5,kphase=6,rout=2.,nz=128):
    psi=potential(H);e0=H.Qt/H.rt;de=e0/M;Npos=int(np.floor(H.psi0/de+.5));Nneg=int(np.ceil(H.psi0/de));ss=(np.arange(-Nneg,Npos)+.5)*de;NE=len(ss);rmax=rout*H.rt;dp=np.sqrt(H.psi0-ss[0])/K
    R=np.r_[0,np.geomspace(.03,rmax,nrad+1),H.rt,rmax];R=np.unique(np.clip(np.r_[R,[brentq(lambda r:float(psi(r))-s,0,rmax,xtol=1e-12) for s in ss if float(psi(rmax))<s<H.psi0]],0,rmax));R=R[np.r_[True,np.diff(R)>1e-10]];NR=len(R)-1
    tab=np.full((NR,NE,2*K+1),-1,int);le=[];lp=[];shell=[];expected=[];ranges=[]
    for j,(a,b) in enumerate(zip(R[:-1],R[1:])):
        ps=float(psi((.5*(a**3+b**3))**(1/3)));v=(b**3-a**3)/3;active=np.flatnonzero(ss<=float(psi(b))+1e-10)
        for i in active:
            E=ps-ss[i];mom=hats(np.sqrt(E),dp,K);sel=np.flatnonzero(mom[:,0]>1e-16);tab[j,i,sel]=np.arange(len(le),len(le)+len(sel));le.extend([-ss[i]]*len(sel));lp.extend((sel-K)*dp);expected.extend(v*de*mom[sel,0]/4);shell.extend([j]*len(sel))
        ranges.append((active,ps,v))
    le=np.array(le);lp=np.array(lp);shell=np.array(shell);expected=np.array(expected);nlocal=len(le);vols=[];oe=[];op=[];oL=[];orbits=[];omass=[];rows=[];cols=[];data=[];dens=[];outside=[];x,w=leggauss(6);xl,wl=leggauss(nL);nopen=0;nbound=0
    for i,s in enumerate(ss):
        if s>=H.psi0:continue
        rs=brentq(lambda r:float(psi(r))-s,0,rmax,xtol=1e-12) if s>float(psi(rmax)) else rmax
        g=lambda r:2*r*r*(float(psi(r))-s);res=minimize_scalar(lambda r:-g(r),bounds=(0,rs),method='bounded',options={'xatol':1e-9});rc=res.x;Lmax=-res.fun
        if g(rs)>Lmax:rc=rs;Lmax=g(rs)
        Lopen=max(0,g(rmax)) if rs==rmax else 0.;Lopen=min(Lopen,Lmax)
        Lcuts=np.unique(np.clip(np.r_[0,Lopen,Lmax,[g(r) for r in R if r<rs]],0,Lmax));Lcuts=Lcuts[np.r_[True,np.diff(Lcuts)>1e-12]]
        for a,b in zip(Lcuts[:-1],Lcuts[1:]):
            isopen=(a+b)/2<Lopen
            if b-a<1e-12:continue
            for z,ww in zip(xl,wl):
                L=(a+b)/2+(b-a)*z/2;WL=(b-a)*ww/2;O=Orbit(psi,s,L,rs,rc,isopen,nz);nopen+=int(isopen);nbound+=int(not isopen)
                ed=mesh(nphase,kphase)
                if isopen:ed=np.r_[0,np.sort(np.mod(ed[1:-1]+.5,1)),1];ed=np.unique(ed)
                if isopen:ed=np.unique(np.r_[ed,.5]);assert len(ed)==nphase+1
                d=np.diff(ed);mid=(ed[:-1]+ed[1:])/2;th=mid[:,None]+d[:,None]*x/2;r,p=O.coords(th);base=len(vols);mass=de*WL*O.T/(8*np.sqrt(2));vols.extend(mass*d);oe.extend([-s]*nphase);oL.extend([L]*nphase);op.extend(p@w/2);orbits.append((O,ed,base));omass.append(mass)
                dens.extend(np.stack([(np.exp(-.5*(r/l)**2)@w/2) for l in [.5,1,2,5,10]],axis=1));frac=np.zeros(nphase)
                for a0,b0 in O.intervals(H.rt,rmax):frac+=np.maximum(0,np.minimum(ed[1:],b0)-np.maximum(ed[:-1],a0))/d
                outside.extend(frac)
                for j,(ra,rb) in enumerate(zip(R[:-1],R[1:])):
                    if tab[j,i].max()<0:continue
                    ps=ranges[j][1];E=ps-s
                    for a0,b0 in O.intervals(ra,rb):
                        ii=np.arange(max(0,np.searchsorted(ed,a0,side='right')-1),min(nphase,np.searchsorted(ed,b0,side='left')));aa=np.maximum(a0,ed[ii]);bb=np.minimum(b0,ed[ii+1]);good=bb>aa;ii=ii[good];aa=aa[good];bb=bb[good];th=(aa[:,None]+bb[:,None])/2+(bb-aa)[:,None]*x/2;r,p=O.coords(th);q=np.sqrt(np.maximum(psi(r)-s,1e-300));val=np.clip(p/q,-1,1)*np.sqrt(E)/dp;ix=np.floor(val).astype(int);f=val-ix;wt=mass*(bb-aa)[:,None]*w/2
                        for k,zwt in [(ix,1-f),(ix+1,f)]:
                            mask=zwt>1e-15;loc=tab[j,i,k[mask]+K];assert np.all(loc>=0);rows.extend(loc);cols.extend(np.broadcast_to((base+ii)[:,None],th.shape)[mask]);data.extend((wt*zwt)[mask])
    Vv=np.array(vols);oe=np.array(oe);op=np.array(op);oL=np.array(oL);dens=np.array(dens);outside=np.array(outside);W=coo_matrix((data,(rows,cols)),shape=(nlocal,len(Vv))).tocsr();local_mass=np.asarray(W.sum(1)).ravel();missing=local_mass<=1e-24
    cerr=float(np.max(abs(np.asarray(W.sum(0)).ravel()/Vv-1)));assert cerr<1e-10,cerr
    samples=qmc.Sobol(4,scramble=True,seed=0).random_base2(m);allids=[];allw=[];lost=0.;total=0.;rho=V**2/(4*np.pi*6.6743e-11*RCONST**2);pref=rho*.125*np.sqrt(2)*V*GYR/H.Ic
    for j,(active,ps,vr) in enumerate(ranges):
        n=len(active)
        for s in range(2*n-1):
            pairs=np.arange(max(0,s-n+1),min(n-1,s//2)+1);other=s-pairs
            for a,b in zip(*np.triu_indices(len(pairs),0)):
                es=active[np.array([pairs[a],other[a],pairs[b],other[b]])];E=ps-ss[es];momidx,v=rounded_distribution(E,dp,K,samples);ids=tab[j,es[:,None],momidx+K];assert np.all(ids>=0)
                wq=pref*vr*de**3*np.sqrt(E).min()/8*(1 if es[0]==es[1] else 2)*(1 if es[2]==es[3] else 2)*(.5 if a==b else 1.)*v/v.sum();keep=~np.any(missing[ids],axis=0);lost+=wq[~keep].sum();total+=wq.sum();allids.append(ids[:,keep].astype(np.int32));allw.append(wq[keep])
    C=Coupling(Vv,oe,W,le,lp,np.concatenate(allids,axis=1),np.concatenate(allw),op);C.shell=shell;C.orbits=orbits;C.omass=np.array(omass);C.dens=dens;C.outside=outside;C.oL=oL;C.e0=e0;C.pref=pref;C.nphase=nphase;C.R=R;C.f0=occ(H,-C.e)
    return C,dict(model=H.model,delta=H.delta,rt=H.rt,e0=e0,de=de,Npos=Npos,Nneg=Nneg,M=M,K=K,nL=nL,nphase=nphase,nrad=nrad,m=m,kphase=kphase,nz=nz,rout=rout,radial_shells=NR,phase_states=len(Vv),local_states=nlocal,events=len(C.weights),closed=nbound,open=nopen,missing_local=int(missing.sum()),lost_event_weight=float(lost/total),overlap_ledger=cerr,volume_L1=float(abs(local_mass-expected).sum()/expected.sum()),pref=pref,initial_mass=float(C.V@C.f0),continuum_mass=H.Ic*H.Qt)

RCONST=R

def run(H,M,K,nL,nphase,nrad,m,kphase,rout,nz,tmax,dt,tilt):
    C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz);f=C.f0.copy();mass=C.V@f;ref=reference(H);ref2=reference(H,32,128);g,b=C.derivative(f,True);e=C.e;ref2*=C.pref/(H.Ic*H.Qt);ref*=C.pref/(H.Ic*H.Qt)
    P=Tilt(C) if tilt else C
    g2,b2=P.derivative(f,True);d.update(reference=ref.tolist(),reference_fine=ref2.tolist(),raw_rate=[float(C.V[e>0]@g[e>0]/mass),float(C.V[e>C.e0]@g[e>C.e0]/mass)],selected_rate=[float(C.V[e>0]@g2[e>0]/mass),float(C.V[e>C.e0]@g2[e>C.e0]/mass)],raw=b,selected=b2,raw_source_L1=float(C.V@abs(g)/mass),correction_L1=float(C.V@abs(g2-g)/mass))
    def stream(f,dt):
        z=f.copy();out=np.zeros(3);minval=0.
        for (O,ed,base),vol in zip(C.orbits,C.omass):
            sl=slice(base,base+nphase);v=f[sl];a=dt/(TIME*O.T)
            if O.open:
                y,l=open_shift(ed,v,a);out+=vol*l*np.array([1,-O.s,float(-O.s>C.e0)])
            else:y=shift(ed,v,a)
            minval=min(minval,float(y.min()));assert y.min()>-1e-9;y=np.maximum(y,0);z[sl]=y
        return z,out,minval
    def rep(f,t,out):
        a=C.stats(f);a.update(t=t,out_count=float(out[0]/mass),out_energy=float(out[1]/mass),out_true_unbound=float(out[2]/mass),cutoff_population=float(C.V[e>0]@f[e>0]/mass),unbound_population=float(C.V[e>C.e0]@f[e>C.e0]/mass),outside_population=float((C.V*C.outside)@f/mass),gaussian=(C.dens.T@(C.V*f)/mass).tolist(),number_ledger=float((C.V@f+out[0])/mass-1),energy_ledger=float(((C.V*e)@f+out[1]-(C.V*e)@C.f0)/mass));return a
    t=0.;out=np.zeros(3);steps=0;reject=0;min_before=0.;d['rows']=[rep(f,t,out)]
    while t<tmax-1e-12:
        h=min(dt,tmax-t)
        while True:
            x,o,mn=stream(f,h/2);g=P.derivative(x);p=x+h*g
            if p.min()>=0:
                y=.5*(x+p+h*P.derivative(p));ds=C.entropy(y)-C.entropy(x) if y.min()>=0 else -1.
                if y.min()>=0 and ds>=-1e-11:break
            h*=.5;reject+=1
            if h<1e-10:raise RuntimeError('collision step below floor')
        f,oo,mn2=stream(y,h/2);out+=o+oo;min_before=min(min_before,mn,mn2);steps+=1;t+=h
    d.update(steps=steps,rejects=reject,dt=dt,tilt=tilt,t=t,min_before_stream_clip=min_before);d['rows'].append(rep(f,t,out));return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=6);p.add_argument('--nphase',type=int,default=64);p.add_argument('--nrad',type=int,default=12);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=6);p.add_argument('--rout',type=float,default=2);p.add_argument('--nz',type=int,default=128);p.add_argument('--t',type=float,default=.02);p.add_argument('--dt',type=float,default=.01);p.add_argument('--tilt',action='store_true');a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt)))
