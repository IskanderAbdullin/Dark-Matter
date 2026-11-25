"""Positive flux-limited equilibrium-subtracted Hamiltonian step.
Conserves number and L including boundary flux. Energy not imposed.
"""
import argparse,json,numpy as np
from lowered_bose_halo_action import Field,build,Halo,TIME,Streamer,Tilt
from lowered_bose_halo_energy_ledger import Energy
class WorkField(Field):
 def __init__(self,H,C,nz=64,mode="work"):
    super().__init__(H,C,nz);self.mode=mode
    if mode=="adjoint":self.energy=Energy(self)
    rr=np.array([o.coords(ed)[0] for o,ed,_ in C.orbits]);rr=np.clip(rr,0,C.rmax);self.cj=np.clip(np.searchsorted(self.r,rr,side='right')-1,0,self.ns-1);bb=self.r[self.cj+1];self.cI=np.divide(1,rr,out=np.zeros_like(rr),where=rr>0)-1/bb;self.cJ=(bb*bb-rr*rr)/2
    self.interior=self.valid_l&self.valid_r;has_l=self.valid_l[self.lo];has_r=self.valid_r[self.hi];wl=self.el*has_l;wr=self.er*has_r;den=wl+wr;self.work_l=np.divide(wl,den,out=np.zeros_like(wl),where=den>0)*C.energy_width;self.work_r=np.divide(wr,den,out=np.zeros_like(wr),where=den>0)*C.energy_width;self.isolated=den==0
    self.max_work_resid=0.;self.max_work_relative=0.;self.work_total=0.;self.work_isolated=0.;self.work_samples=0
 def center_potential(self,q):
    a,b=self.shell_a,self.shell_b;c1=np.diff(q)/self.shell_vol;c0=q[:-1]-c1*a**3;ia=np.divide(1,a,out=np.zeros_like(a),where=a>0)-1/b;I=c0*ia+c1*(b*b-a*a)/2;suf=np.r_[np.cumsum(I[::-1])[::-1],0.];j=self.cj
    return -q[-1]/self.r[-1]-suf[j+1]-c0[j]*self.cI-c1[j]*self.cJ
 def current_flux(self,z,q,E0):
    if self.mode=="adjoint":
     f=z.reshape(self.no,self.n);phi=(self.energy.hamiltonian(z+self.C.f0)-self.C.e).reshape(self.no,self.n)
     left=np.roll(f,1,axis=1);right=np.roll(f,-1,axis=1);a=(f-left)/(self.tc-self.tl);b=(right-f)/(self.tr-self.tc);sl=np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.);full=f+self.C.f0.reshape(self.no,self.n);sl=np.clip(sl,-2*full/self.d,2*full/self.d);sl[self.open,0]=sl[self.open,-1]=0.
     face=f+sl*self.d/2;dphi=np.roll(phi,-1,axis=1)-phi;dphi[self.open,-1]=0.;p=face*dphi;prev=np.roll(p,1,axis=1);prev[self.open,0]=0.;B=-self.w[:,None]*(p+prev)/2
    else:
     dh=np.diff(self.center_potential(q),axis=1);f=z.reshape(self.no,self.n);B=-self.w[:,None]*f*dh
    target=float(np.sum(self.C.energy_width[:,None]*B));l=self.lf[self.interior];r=self.rf[self.interior];gap=self.e[r]-self.e[l];E=E0.copy();E[self.interior]=(self.work_r[l,None]*B[l]+self.work_l[r,None]*B[r])/gap[:,None]
    moment=float(np.sum(gap[:,None]*E[self.interior]));iso=float(np.sum(self.C.energy_width[self.isolated,None]*B[self.isolated]));self.max_work_resid=max(self.max_work_resid,abs(moment+iso-target));self.max_work_relative=max(self.max_work_relative,abs(moment+iso-target)/max(abs(target),1e-20));self.work_total+=abs(target);self.work_isolated+=abs(iso);self.work_samples+=1;return E

 def euler(self,f,coll,dt):
    V=self.C.V; q=self.enclosed(f);E0,T0,*_=self.faces(f,q,1);E1,T1,*_=self.faces(f-self.C.f0,q,2,True)
    rev=(f-self.C.f0).reshape(self.no,self.n)[:,::-1].ravel();Er,Tr,*_=self.faces(rev,q,2,True);E1=(E1-Er[:,::-1])/2;Ts=(T1+Tr[:,::-1])/2;Ts[self.open,0]=T0[self.open,0];Ts[self.open,-1]=T0[self.open,-1];T1=Ts;E1[~(self.valid_l&self.valid_r)]=E0[~(self.valid_l&self.valid_r)]
    if self.mode!="parity":E1=self.current_flux(f-self.C.f0,q,E0)
    eq=(V*(self.F.G@q)).reshape(self.no,self.n)
    T1[:,1:]-=np.cumsum(eq,axis=1);T1[~self.open,-1]=T1[~self.open,0]
    g0,o0,*_=self.counts(E0,T0);low=f+dt*(coll+g0.ravel()/V)
    if low.min()<-1e-14:return None,None,None
    clipped=float(V@np.maximum(-low,0));low=np.maximum(low,0)
    dE=E1-E0;dT=T1-T0
    outgoing=np.maximum(dE[self.hi],0)-np.minimum(dE[self.lo],0)+np.maximum(dT[:,1:],0)-np.minimum(dT[:,:-1],0)
    available=(V*low).reshape(self.no,self.n);ratio=np.minimum(1,np.divide(available,dt*outgoing,out=np.ones_like(available),where=outgoing>0))
    a=np.ones_like(dE);b=np.ones_like(dE);a[self.valid_l]=ratio[self.lf[self.valid_l]];b[self.valid_r]=ratio[self.rf[self.valid_r]];se=np.where(dE>=0,a,b)
    a=np.hstack((ratio[:,-1:],ratio));b=np.hstack((ratio,ratio[:,:1]));a[self.open,0]=1;b[self.open,-1]=1;st=np.where(dT>=0,a,b)
    st[~self.open,-1]=st[~self.open,0]
    ge,oe,*_=self.counts(se*dE,st*dT);result=low+dt*ge.ravel()/V
    if result.min()<-1e-13:raise RuntimeError(('limiter positivity',result.min()))
    clipped+=float(V@np.maximum(-result,0));result=np.maximum(result,0)
    raw=float(abs(dE).sum()+abs(dT).sum());limited=float((abs(dE)*(1-se)).sum()+(abs(dT)*(1-st)).sum())
    return result,dt*(o0+oe),dict(raw=dt*raw,limited=dt*limited,clip=clipped,eq_ledger=float(abs(eq.sum(1)).max()))

def run(H,M=1,K=4,nL=1,nphase=32,nrad=10,m=5,kphase=4,rout=2.,nz=64,tmax=.1,dt=.005,tilt=False,mode="work"):
 C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz,True);P=Tilt(C) if tilt else C;F=WorkField(H,C,nz,mode);S=Streamer(C,TIME);V=C.V;e=C.e;f0=C.f0.copy();f=f0.copy();fr=f0.copy();mass=V@f0;out=np.zeros(2);outr=np.zeros(2);steps=rej=0;t=0.;diag=dict(raw=0.,limited=0.,clip=0.,eq_ledger=0.)
 def step(f,h,gr):
    x,o,_=S.apply(f,h/2)
    if gr:
        p,oa,da=F.euler(x,P.derivative(x),h)
        if p is None:return None,None,None
        y,ob,db=F.euler(p,P.derivative(p),h)
        if y is None:return None,None,None
        y=(x+y)/2;dd={k:(max(da[k],db[k]) if k=='eq_ledger' else (da[k]+db[k])/2) for k in da}
    else:
        p=x+h*P.derivative(x)
        if p.min()<-1e-14:return None,None,None
        p=np.maximum(p,0);y=(x+p+h*P.derivative(p))/2
        if y.min()<-1e-14:return None,None,None
        y=np.maximum(y,0);oa=ob=np.zeros(3);dd=None
    z,oo,_=S.apply(y,h/2);return z,o[:2]+oo[:2]+(oa[:2]+ob[:2])/2,dd
 def report(f,t,out):
    field=F.F.report(f);E=float(((V*e)@(f-f0)+out[1])/mass);return dict(t=t,Nledger=float((V@f+out[0])/mass-1),reference_energy_change=E,field_energy=field['perturbation_field_energy'],total_energy_error=E+field['perturbation_field_energy'],min=float(f.min()),max_force_rel=field['max_relative_force'],gaussian_change=(C.dens.T@(V*(f-f0))/mass).tolist(),cutoff=float(V[e>0]@f[e>0]/mass),unbound=float(V[e>C.e0]@f[e>C.e0]/mass),outside=float((V*C.outside)@f/mass),out_count=float(out[0]/mass))
 rows=[dict(full=report(f,0,out),frozen=report(fr,0,outr))]
 while t<tmax-1e-12:
    h=min(dt,tmax-t)
    while True:
        a,b,dd=step(f,h,True);c,ddd,_=step(fr,h,False)
        if a is not None and c is not None:break
        h*=.5;rej+=1
        if h<1e-8:raise RuntimeError('positivity failure')
    f,fr=a,c;out+=b;outr+=ddd;t+=h;steps+=1
    for k in diag:diag[k]=max(diag[k],dd[k]) if k=='eq_ledger' else diag[k]+dd[k]
 rows.append(dict(full=report(f,t,out),frozen=report(fr,t,outr)));d.update(mode=mode,work_diagnostic=dict(max_residual=F.max_work_resid,max_relative=F.max_work_relative,isolated_relative=F.work_isolated/max(F.work_total,1e-300),samples=F.work_samples),steps=steps,rejections=rej,rows=rows,limiter=diag,relative_limited_flux=diag['limited']/max(diag['raw'],1e-300),dt=dt,tmax=tmax);return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=1);p.add_argument('--nphase',type=int,default=32);p.add_argument('--nrad',type=int,default=10);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=4);p.add_argument('--rout',type=float,default=2);p.add_argument('--nz',type=int,default=64);p.add_argument('--t',type=float,default=.1);p.add_argument('--dt',type=float,default=.005);p.add_argument('--tilt',action='store_true');p.add_argument('--mode',choices=['parity','work','adjoint'],default='work');a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt,a.mode)))
