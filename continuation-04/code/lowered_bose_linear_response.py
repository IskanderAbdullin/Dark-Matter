"""Causal spherical linearized Vlasov--Poisson response to constant C22[f0].
The collision operator is not linearized beyond its initial forcing. Each
field-response kernel column is an independently streamed odd-in-v_r source.
"""
import argparse,json,math
import numpy as np
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes
from lowered_bose_signed_response import source
from lowered_bose_exterior_source import TIME

def equilibrium_kernel(H,nr,nq,na):
    x,w=nodes(nr);r=H.rt*x*x;dr=w*2*H.rt*x;ps=H.psi(r)
    z,wq=nodes(nq);mu,wa=nodes(na);mu=2*mu-1
    gap=H.eta-ps;thmax=np.arctan(np.sqrt(ps/gap));th=thmax[:,None]*z
    q=np.sqrt(gap[:,None])*np.tan(th);dq=np.sqrt(gap[:,None])/np.cos(th)**2*thmax[:,None]*wq
    s=ps[:,None]-q*q;b=1/np.expm1(H.eta-s);Fp=b*(1+b)
    if H.model=='quadratic':Fp=np.where(s<1e-5,H.d*(1+2*H.c)*s+H.d*(1+6*H.c+6*H.c**2)*s*s/2,Fp-H.d)
    F=H.F(s);sus=(dq*q*q*Fp).sum(axis=1);other=.5*(dq*F).sum(axis=1)
    pp=ps[:,None];qq=np.sqrt(pp)*H.u;vv=gap[:,None]+qq*qq;ss=pp-qq*qq
    reg=np.where(vv<1e-3,-.5+vv/12-vv**3/720,1/np.expm1(vv)-1/vv)-H.c
    if H.model=='quadratic':reg-=H.d*ss
    exact=.5*(np.arctan(np.sqrt(ps/gap))/np.sqrt(gap)+np.sqrt(ps)*np.sum(H.w*reg,axis=1))
    exact=np.where(ps<.2,.5*np.sqrt(ps)*np.sum(H.w*H.F(ss),axis=1),exact)
    vr=np.sqrt(2)*q[:,:,None]*mu;vt=np.sqrt(2)*q[:,:,None]*np.sqrt(1-mu*mu)
    W=-dr[:,None,None]*(dq*q*q*Fp)[:,:,None]*vr*wa/(H.Ic*TIME)
    group=np.broadcast_to(np.arange(nr)[:,None,None],vr.shape).ravel()
    arr=np.stack(np.broadcast_arrays(np.broadcast_to(r[:,None,None],vr.shape),vr,vt,W),axis=-1).reshape(-1,4)
    rel=abs(sus/exact-1);weighted=np.sum(dr*r*r*abs(sus-exact))/np.sum(dr*r*r*exact)
    return r,arr,group,dict(sus_max=float(rel.max()),sus_weighted=float(weighted),sus_ibp_weighted=float(np.sum(dr*r*r*abs(sus-other))/np.sum(dr*r*r*exact)))

def smooth_enclosed(r,w,group,mid,ng):
    nr=len(mid)-1;idx=np.clip(np.searchsorted(mid,r,side='right')-1,0,nr);flat=group*(nr+1)+idx
    counts=np.bincount(flat,weights=w,minlength=ng*(nr+1)).reshape(ng,nr+1)[:,:nr]
    j=np.minimum(idx,nr-1);frac=np.where(idx<nr,np.clip((mid[j+1]-r)/(mid[j+1]-mid[j]),0,1),0)
    part=np.bincount(flat,weights=w*frac,minlength=ng*(nr+1)).reshape(ng,nr+1)[:,:nr]
    return (np.cumsum(counts,axis=1)-counts+part).T

def volume_windows(rgrid,rt):
    y=rgrid**3;prev=np.r_[0.,y[:-1]];nxt=np.r_[y[1:],rt**3]
    dy=.5*np.minimum(y-prev,nxt-y)
    return y-dy,y+dy

def volume_enclosed(r,w,group,window,ng):
    lo,hi=window;nr=len(lo);y=r**3;idx=np.searchsorted(hi,y,side='right');flat=group*(nr+1)+idx
    counts=np.bincount(flat,weights=w,minlength=ng*(nr+1)).reshape(ng,nr+1)[:,:nr]
    j=np.minimum(idx,nr-1);frac=np.where(idx<nr,np.clip((hi[j]-y)/(hi[j]-lo[j]),0,1),0)
    part=np.bincount(flat,weights=w*frac,minlength=ng*(nr+1)).reshape(ng,nr+1)[:,:nr]
    return (np.cumsum(counts,axis=1)-counts+part).T

def solve(H,rgrid,K,born,h,scales,radii,scenarios):
    nr=len(rgrid);steps=len(born)-1;ns=len(scenarios);d=np.zeros((ns,steps+1,born.shape[1]));fields=np.zeros((ns,steps+1,nr));anchors=[]
    for rc,p in scenarios:
        j=int(np.searchsorted(rgrid,rc));anchors.append(j)
    marks={round(t/h):t for t in [.01,.05,.1,.2,.5,1.] if t<=steps*h+1e-12 and round(t/h)>0};rows=[[] for _ in scenarios]
    base=np.interp(rgrid,np.r_[0.,rgrid,H.rt],np.r_[0.,H.sol.sol(rgrid)[1],H.Qt])/rgrid**2
    for i in range(1,steps+1):
        d[:,i]=born[i]+h*np.einsum('tij,stj->si',K[1:i+1],fields[:,i-1::-1],optimize=False)
        fields[:,i]=d[:,i,:nr]
        for j,((rc,p),a) in enumerate(zip(scenarios,anchors)):
            if a:fields[j,i,:a]=fields[j,i,a]*(rgrid[:a]/rgrid[a])**p
        if i in marks or i==steps:
            sl=slice(nr,nr+len(scales));en=slice(nr+len(scales),nr+len(scales)+len(radii));mask=rgrid>=.1
            for j in range(ns):
                x=d[j,i];f=H.Qt*fields[j,i]/rgrid**2
                rows[j].append(dict(time_Gyr=i*h,born_smooth=born[i,sl].tolist(),response_smooth=x[sl].tolist(),born_enclosed=born[i,en].tolist(),response_enclosed=x[en].tolist(),number_ledger=float(x[-2]),energy_ledger=float(x[-1]),max_enclosed_grid=float(np.max(abs(x[:nr]))),max_raw_force_grid=float(np.max(abs(H.Qt*x[:nr]/rgrid**2))),max_effective_force_grid=float(np.max(abs(f))),max_relative_force_outside_01=float(np.max(abs(f[mask]/base[mask])))))
    return [dict(core_requested=rc,power=p,anchor_radius=float(rgrid[a]) if a else 0.,rows=row) for (rc,p),a,row in zip(scenarios,anchors,rows)]

def run(H,n=12,nr=32,nq=32,na=10,dt=.002,T=.5,window='volume',multi=True):
    src,source_err=source(H,n,nr,na);rgrid,ker,kg,kcheck=equilibrium_kernel(H,nr,nq,na)
    a=np.concatenate((src,ker));x,u,v,w=a.T.copy();y=np.zeros_like(x);nc=len(src);nk=len(ker);wc=w[:nc];wk=w[nc:];cg=np.zeros(nc,dtype=int)
    rt=H.rt;e0=H.Qt/rt;rr=np.r_[0.,rt*np.geomspace(1e-7,1,6000)];mass=H.sol.sol(np.maximum(rr,H.eps))[1];coeff=np.where(rr>H.eps,mass/np.maximum(rr,1e-30)**3,1/3);U=e0+H.psi(rr)
    def force(r):return np.where(r<=rt,np.interp(r,rr,coeff),H.Qt/np.maximum(r,1e-30)**3)
    def potential(r):return np.where(r<=rt,np.interp(r,rr,U),H.Qt/np.maximum(r,1e-30))
    def energy(r):return .5*(u*u+v*v)-potential(r)
    E0=energy(x);mid=np.r_[0.,(rgrid[1:]+rgrid[:-1])/2,rt]
    win=volume_windows(rgrid,rt);enclose=(lambda r,w,g,ng:volume_enclosed(r,w,g,win,ng)) if window=="volume" else (lambda r,w,g,ng:smooth_enclosed(r,w,g,mid,ng))
    scales=np.array([.5,1.,2.,5.,10.]);radii=np.array([1.,2.,5.,10.,rt,2*rt,5*rt]);nobs=nr+len(scales)+len(radii)+2
    def obs(r):
        cr=r[:nc];kr=r[nc:]
        cs=enclose(cr,wc,cg,1)[:,0];ks=enclose(kr,wk,kg,nr)
        cs=np.r_[cs,np.exp(-.5*(cr[None,:]/scales[:,None])**2)@wc,[wc[cr<rad].sum() for rad in radii],wc.sum(),np.dot(wc,energy(r)[:nc])]
        kk=[np.bincount(kg,weights=wk*np.exp(-.5*(kr/l)**2),minlength=nr) for l in scales]
        kk+=[np.bincount(kg,weights=wk*(kr<rad),minlength=nr) for rad in radii]
        kk+=[np.bincount(kg,weights=wk,minlength=nr),np.bincount(kg,weights=wk*energy(r)[nc:],minlength=nr)]
        return cs,np.vstack((ks,kk))
    r=x.copy();prev,k0=obs(r);initial_cs=float(np.max(abs(prev)));initial_k=float(np.max(abs(k0)));steps=round(T/dt);h=T/steps;ht=h/TIME
    K=np.zeros((steps+1,nobs,nr));born=np.zeros((steps+1,nobs));f=force(r);ax=-f*x;ay=-f*y
    for i in range(1,steps+1):
        u+=.5*ht*ax;v+=.5*ht*ay;x+=ht*u;y+=ht*v;r=np.hypot(x,y);f=force(r);ax=-f*x;ay=-f*y;u+=.5*ht*ax;v+=.5*ht*ay
        cur,K[i]=obs(r);born[i]=born[i-1]+.5*h*(prev+cur);prev=cur
    scenarios=[(0.,3),(.05,3),(.1,3),(.2,3),(.1,2)] if multi else [(0.,3)]
    results=solve(H,rgrid,K,born,h,scales,radii,scenarios)
    return dict(model=H.model,delta=H.delta,n=n,nr=nr,nq=nq,na=na,dt=h,window=window,source_tracers=nc,kernel_tracers=nk,initial_source_obs=initial_cs,initial_kernel_obs=initial_k,kernel_checks=kcheck,source_ledger_max=source_err.tolist(),kernel_number_max=float(np.max(abs(K[:,-2]))),kernel_energy_max=float(np.max(abs(K[:,-1]))),scenarios=results)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=0);p.add_argument('--n',type=int,default=12);p.add_argument('--nr',type=int,default=32);p.add_argument('--nq',type=int,default=32);p.add_argument('--na',type=int,default=10);p.add_argument('--dt',type=float,default=.002);p.add_argument('--T',type=float,default=.5);p.add_argument('--window',choices=['linear','volume'],default='volume');p.add_argument('--single',action='store_true');a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model,delta=a.delta),a.n,a.nr,a.nq,a.na,a.dt,a.T,a.window,not a.single)))
