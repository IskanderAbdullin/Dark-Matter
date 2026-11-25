"""Regular radial weak-Poisson Galerkin causal VP response (first order C22)."""
import argparse,json,numpy as np
from scipy.linalg import cho_factor,cho_solve
from lowered_antonov_hessian import Halo
from lowered_bose_signed_source import nodes
from lowered_bose_signed_response import source
from lowered_bose_exterior_source import TIME

def radial_basis(H,nf,r0,ng):
    grid=np.r_[0,np.geomspace(r0,H.rt,nf)];a,b=grid[:-1],grid[1:];den=b*b-a*a
    val=4*(b**5-a**5)/(5*den**2);A=np.zeros((nf,nf))
    for i,c in enumerate(val):
        A[i,i]+=c
        if i+1<nf:A[i+1,i+1]+=c;A[i,i+1]-=c;A[i+1,i]-=c
    x,w=nodes(ng);r=(a[:,None]+(b-a)[:,None]*x).ravel();dr=((b-a)[:,None]*w).ravel();k=np.repeat(np.arange(nf),ng)
    D=np.zeros((len(r),nf+1));j=np.arange(len(r));D[j,k]=2*r**3/(H.Qt*den[k]);D[j,k+1]=-D[j,k]
    return grid,A,r,dr,D[:,:nf]

def kernel(H,r,dr,nq,na):
    ps=H.psi(r);z,wq=nodes(nq);mu,wa=nodes(na);mu=2*mu-1
    gap=H.eta-ps;thmax=np.arctan(np.sqrt(ps/gap));th=thmax[:,None]*z
    q=np.sqrt(gap[:,None])*np.tan(th);dq=np.sqrt(gap[:,None])/np.cos(th)**2*thmax[:,None]*wq
    s=ps[:,None]-q*q;b=1/np.expm1(H.eta-s);fp=b*(1+b)
    if H.model=='quadratic':fp-=H.d
    vr=np.sqrt(2)*q[:,:,None]*mu;vt=np.sqrt(2)*q[:,:,None]*np.sqrt(1-mu*mu)
    W=-dr[:,None,None]*(dq*q*q*fp)[:,:,None]*vr*wa/(H.Ic*TIME)
    group=np.broadcast_to(np.arange(len(r))[:,None,None],vr.shape).ravel()
    arr=np.stack(np.broadcast_arrays(np.broadcast_to(r[:,None,None],vr.shape),vr,vt,W),axis=-1).reshape(-1,4)
    return arr,group

def weak_mass(r,w,g,ng,grid):
    nf=len(grid)-1;k=np.clip(np.searchsorted(grid,r,side='right')-1,0,nf-1);t=np.clip((r*r-grid[k]**2)/(grid[k+1]**2-grid[k]**2),0,1)
    f=g*(nf+1)+k
    a=np.bincount(f,weights=w*(1-t),minlength=ng*(nf+1));b=np.bincount(f+1,weights=w*t,minlength=ng*(nf+1));return (a+b).reshape(ng,nf+1)[:,:nf].T

def run(H,n=12,nr=32,nq=24,na=10,dt=.002,T=.5,nf=24,r0=.25,ng=4):
    grid,A,rg,dw,D=radial_basis(H,nf,r0,ng);fac=cho_factor(A);src,source_err=source(H,n,nr,na);ker,kg=kernel(H,rg,dw,nq,na);nk=len(rg)
    a=np.concatenate((src,ker));x,u,v,w=a.T.copy();y=np.zeros_like(x);nc=len(src);wc=w[:nc];wk=w[nc:];cg=np.zeros(nc,dtype=int)
    rt=H.rt;e0=H.Qt/rt;rr=np.r_[0.,rt*np.geomspace(1e-7,1,6000)];mass=H.sol.sol(np.maximum(rr,H.eps))[1];coeff=np.where(rr>H.eps,mass/np.maximum(rr,1e-30)**3,1/3);U=e0+H.psi(rr)
    def force(r):return np.where(r<=rt,np.interp(r,rr,coeff),H.Qt/np.maximum(r,1e-30)**3)
    def potential(r):return np.where(r<=rt,np.interp(r,rr,U),H.Qt/np.maximum(r,1e-30))
    def energy(r):return .5*(u*u+v*v)-potential(r)
    scales=np.array([.5,1.,2.,5.,10.]);radii=np.array([1.,2.,5.,10.,rt,2*rt,5*rt]);nobs=nf+len(scales)+len(radii)+2
    def obs(r):
        cr=r[:nc];kr=r[nc:]
        cs=weak_mass(cr,wc,cg,1,grid)[:,0];ks=weak_mass(kr,wk,kg,nk,grid)
        cs=np.r_[cs,np.exp(-.5*(cr[None,:]/scales[:,None])**2)@wc,[wc[cr<rad].sum() for rad in radii],wc.sum(),np.dot(wc,energy(r)[:nc])]
        kk=[np.bincount(kg,weights=wk*np.exp(-.5*(kr/l)**2),minlength=nk) for l in scales]
        kk+=[np.bincount(kg,weights=wk*(kr<rad),minlength=nk) for rad in radii]
        kk+=[np.bincount(kg,weights=wk,minlength=nk),np.bincount(kg,weights=wk*energy(r)[nc:],minlength=nk)]
        return cs,np.vstack((ks,kk))@D
    r=x.copy();prev,k0=obs(r);steps=round(T/dt);h=T/steps;ht=h/TIME;K=np.zeros((steps+1,nobs,nf));born=np.zeros((steps+1,nobs));f=force(r);ax=-f*x;ay=-f*y
    for i in range(1,steps+1):
        u+=.5*ht*ax;v+=.5*ht*ay;x+=ht*u;y+=ht*v;r=np.hypot(x,y);f=force(r);ax=-f*x;ay=-f*y;u+=.5*ht*ax;v+=.5*ht*ay
        cur,K[i]=obs(r);born[i]=born[i-1]+.5*h*(prev+cur);prev=cur
    d=np.zeros_like(born);field=np.zeros((steps+1,nf));rows=[];marks={round(t/h):t for t in [.01,.05,.1,.2,.5,1.] if t<=steps*h+1e-12 and round(t/h)>0}
    base=H.sol.sol(rg)[1]/rg**2;mask=rg>=.1
    for i in range(1,steps+1):
        d[i]=born[i]+h*np.einsum('tij,tj->i',K[1:i+1],field[i-1::-1],optimize=False);field[i]=cho_solve(fac,H.Qt*d[i,:nf]);M=D@field[i];f=H.Qt*M/rg**2
        if i in marks or i==steps:
            sl=slice(nf,nf+len(scales));en=slice(nf+len(scales),nf+len(scales)+len(radii));rows.append(dict(time_Gyr=i*h,born_smooth=born[i,sl].tolist(),response_smooth=d[i,sl].tolist(),response_enclosed=d[i,en].tolist(),number_ledger=float(d[i,-2]),energy_ledger=float(d[i,-1]),max_force=float(np.max(abs(f))),max_rel_force_outside_01=float(np.max(abs(f[mask]/base[mask]))),central_force_slope=float(2*(field[i,0]-field[i,1])/grid[1]**2),poisson_residual=float(np.max(abs(A@field[i]-H.Qt*d[i,:nf])))))
    return dict(model=H.model,delta=H.delta,n=n,nr=nr,nq=nq,na=na,dt=h,nf=nf,r0=r0,ng=ng,source_tracers=nc,kernel_tracers=len(ker),initial_kernel_max=float(np.max(abs(k0))),source_ledger_max=source_err.tolist(),kernel_number_max=float(np.max(abs(K[:,-2]))),rows=rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=0);p.add_argument('--n',type=int,default=12);p.add_argument('--nr',type=int,default=32);p.add_argument('--nq',type=int,default=24);p.add_argument('--na',type=int,default=10);p.add_argument('--dt',type=float,default=.002);p.add_argument('--T',type=float,default=.5);p.add_argument('--nf',type=int,default=24);p.add_argument('--r0',type=float,default=.25);p.add_argument('--ng',type=int,default=4);a=p.parse_args();print(json.dumps(run(Halo(model=a.model,delta=a.delta),a.n,a.nr,a.nq,a.na,a.dt,a.T,a.nf,a.r0,a.ng)))
