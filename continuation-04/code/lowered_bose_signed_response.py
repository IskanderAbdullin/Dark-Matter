"""First Born/Duhamel response: stream the entire signed initial C22 source in
frozen potential. No linearized gravity or repeated C22; not self-consistent evolution.
"""
import argparse,json,math
import numpy as np
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes,local
from lowered_bose_exterior_source import R,V,GYR,TIME

def source(H,n,nr,na):
    e0=H.Qt/H.rt;eta=H.eta
    edges=np.unique(np.r_[np.linspace(-eta,-e0,13),np.linspace(-e0,0,13),np.linspace(0,eta-.2,40),eta-np.geomspace(.2,.001,12),eta]);edges[0]-=1e-12;edges[-1]+=1e-12
    x,w=nodes(nr);r=H.rt*x*x;ww=w*2*H.rt*x*r*r;ps=H.psi(r)
    rho=V**2/(4*math.pi*6.6743e-11*R**2);pref=.5*rho*.125*math.sqrt(2)*V*GYR/(H.Ic**2*H.Qt)
    mu,wm=nodes(na);mu=2*mu-1;rows=[];errs=[]
    for rr,p,wgt in zip(r,ps,ww):
        c,m,aw=local(H,p,edges,n)
        lo=edges[:-1];hi=np.minimum(edges[1:],p);mask=(hi>lo)&((c!=0)|(m!=0));lo=lo[mask];hi=hi[mask];c=c[mask];m=m[mask]
        wl=(hi*c-m)/(hi-lo);wh=c-wl
        ss=np.r_[lo,hi];weights=np.r_[wl,wh]*wgt*pref
        errs.append([abs(weights.sum()),abs(np.dot(weights,ss))])
        v=np.sqrt(2*np.maximum(0,p-ss));vr=v[:,None]*mu;vt=v[:,None]*np.sqrt(1-mu*mu);weight=weights[:,None]*wm
        rows.append(np.stack(np.broadcast_arrays(np.full_like(vr,rr),vr,vt,weight),axis=-1).reshape(-1,4))
    a=np.concatenate(rows);return a,np.max(errs,axis=0)

def run(H,n=12,nr=32,na=6,dt=.001,T=1.):
    a,source_err=source(H,n,nr,na);x,u,v,w=a.T.copy();y=np.zeros_like(x)
    rt=H.rt;e0=H.Qt/rt
    rr=np.r_[0.,rt*np.geomspace(1e-7,1,6000)];mass=H.sol.sol(np.maximum(rr,H.eps))[1];coeff=np.where(rr>H.eps,mass/np.maximum(rr,1e-30)**3,1/3);U=e0+H.psi(rr)
    def force(r):return np.where(r<=rt,np.interp(r,rr,coeff),H.Qt/np.maximum(r,1e-30)**3)
    def potential(r):return np.where(r<=rt,np.interp(r,rr,U),H.Qt/np.maximum(r,1e-30))
    def energy(r):return .5*(u*u+v*v)-potential(r)
    E0=energy(x);absw=np.sum(abs(w));wE0=np.dot(w,E0)
    # smooth mass probes and sharp enclosed radii; total response exactly zero.
    scales=np.array([.5,1.,2.,5.,10.]);radii=np.array([1.,2.,5.,10.,rt,2*rt,5*rt]);
    def obs(r):
        return np.r_[np.exp(-.5*(r[None,:]/scales[:,None])**2)@w,np.array([w[r<rad].sum() for rad in radii])]
    r=x.copy();prev=obs(r);integ=np.zeros_like(prev);initial=prev.copy();rows=[];steps=round(T/dt);h=T/steps;ht=h/TIME;f=force(r);ax=-f*x;ay=-f*y
    marks={round(t/h):t for t in [.01,.05,.1,.2,.5,1.,2.,5.] if t<=T+1e-12 and round(t/h)>0}
    for k in range(1,steps+1):
        u+=.5*ht*ax;v+=.5*ht*ay;x+=ht*u;y+=ht*v;r=np.hypot(x,y);f=force(r);ax=-f*x;ay=-f*y;u+=.5*ht*ax;v+=.5*ht*ay
        cur=obs(r);integ+=.5*h*(prev+cur);prev=cur
        if k in marks or k==steps:
            E=energy(r);rows.append(dict(time_Gyr=k*h,smooth=integ[:len(scales)].tolist(),enclosed=integ[len(scales):].tolist(),number_ledger=k*h*w.sum(),weighted_energy_initial=wE0,weighted_energy_drift=np.dot(w,E-E0),mean_abs_energy_drift=np.dot(abs(w),abs(E-E0))/absw))
    return dict(model=H.model,n=n,nr=nr,na=na,dt=h,trajectories=len(w),abs_source_weight=absw,source_ledger_max=source_err.tolist(),scales=scales.tolist(),radii=radii.tolist(),initial_obs=initial.tolist(),rows=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--n',type=int,default=12);p.add_argument('--nr',type=int,default=32);p.add_argument('--na',type=int,default=6);p.add_argument('--dt',type=float,default=.001);p.add_argument('--T',type=float,default=1.);a=p.parse_args()
    print(json.dumps(run(Halo(model=a.model),a.n,a.nr,a.na,a.dt,a.T)))
