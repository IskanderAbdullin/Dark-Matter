"""Diagnostic adjacent-energy work repair, not a derived continuum force."""
import argparse,json,numpy as np
from scipy.optimize import brentq
from lowered_bose_halo_energy_ledger import Energy,LimitedField,build,Halo,TIME,Streamer,Tilt

class Repair:
 def __init__(self,F,E):self.F=F;self.E=E;self.V=E.V;self.inner=F.valid_l&F.valid_r
 def apply(self,start,f,out_energy,dt):
    E,F,V=self.E,self.F,self.V;err=E.difference(start,f,out_energy).sum()*E.mass
    if abs(err)<1e-15:return f,dict(error=err/E.mass,after=err/E.mass,alpha=0.,l1=0.,relative_l1=0.,flux_l1=0.,max_frac=0.,potential_ratio=0.)
    q=F.enclosed(f);EF=F.faces(f-F.C.f0,q,2,True)[0];EF[~self.inner]=0.
    flux=np.maximum(EF,0) if err<0 else np.minimum(EF,0)
    outgoing=np.maximum(flux[F.hi],0)-np.minimum(flux[F.lo],0);available=(V*f).reshape(F.no,F.n);cap=np.divide(available,outgoing,out=np.full_like(available,np.inf),where=outgoing>0)
    def transfer(alpha):
        fac=np.minimum(alpha,cap);a=np.zeros_like(flux);b=np.zeros_like(flux);a[F.valid_l]=fac[F.lf[F.valid_l]];b[F.valid_r]=fac[F.rf[F.valid_r]];used=flux*np.where(flux>=0,a,b);dn=-(used[F.hi]-used[F.lo]).ravel();return dn,used
    def fn(alpha):
        dn,_=transfer(alpha);df=dn/V;dm=E.W@df;return err+(V*E.hamiltonian(f))@df-.5*dm@E.A@dm
    bound=dt
    while fn(bound)*err>0 and bound<dt*1e5:bound*=2
    if fn(bound)*err>0:raise RuntimeError(('repair capacity insufficient',err,bound,fn(bound)))
    alpha=brentq(fn,0,bound,xtol=1e-16,rtol=1e-14);dn,used=transfer(alpha);g=f+dn/V
    if g.min()<-1e-14:raise RuntimeError('repair positivity failure')
    g=np.maximum(g,0);after=E.difference(start,g,out_energy).sum();raw_l1=float(abs(V*(f-start)).sum());out_used=np.maximum(used[F.hi],0)-np.minimum(used[F.lo],0);frac=np.divide(out_used,available,out=np.zeros_like(available),where=available>0);potential=float((E.V*(E.hamiltonian(f)-E.e))@(dn/V));linear=float(E.e@dn)
    return g,dict(error=err/E.mass,after=after,alpha=float(alpha/dt),l1=float(abs(dn).sum()/E.mass),relative_l1=float(abs(dn).sum()/max(raw_l1,1e-300)),flux_l1=float(dt*abs(flux).sum()/E.mass),max_frac=float(frac.max()),potential_ratio=float(potential/max(abs(linear),1e-300)))

def run(H,M=1,K=4,nL=1,nphase=32,nrad=10,m=5,kphase=4,rout=2.,nz=64,tmax=.1,dt=.005,tilt=False):
 C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz,True);P=Tilt(C) if tilt else C;F=LimitedField(H,C,nz);S=Streamer(C,TIME);E=Energy(F);R=Repair(F,E);f=C.f0.copy();fr=f.copy();out=np.zeros(2);outr=np.zeros(2);t=0.;steps=rej=0;history=[]
 def step(f,h):
    x,o,_=S.apply(f,h/2);p,oa,da=F.euler(x,P.derivative(x),h)
    if p is None:return None,None
    y,ob,db=F.euler(p,P.derivative(p),h)
    if y is None:return None,None
    y=(x+y)/2;z,oo,_=S.apply(y,h/2);return z,o[:2]+oo[:2]+(oa[:2]+ob[:2])/2
 while t<tmax-1e-12:
    h=min(dt,tmax-t)
    while True:
        a,b=step(f,h);c,dd=step(fr,h)
        if a is not None and c is not None:break
        h*=.5;rej+=1
        if h<1e-8:raise RuntimeError('positivity failure')
    a,row=R.apply(f,a,b[1],h);row['t']=t+h;history.append(row);f,fr=a,c;out+=b;outr+=dd;t+=h;steps+=1
 def report(f,out):
    ed=E.difference(C.f0,f,out[1]);return dict(Nledger=float((E.V@f+out[0])/E.mass-1),Eref=float(ed[0]),field=float(ed[1]),energy_residual=float(ed.sum()),gaussian=(C.dens.T@(E.V*(f-C.f0))/E.mass).tolist(),cutoff=float(E.V[E.e>0]@f[E.e>0]/E.mass),unbound=float(E.V[E.e>C.e0]@f[E.e>C.e0]/E.mass),out_count=float(out[0]/E.mass),min=float(f.min()),Ltotal=float((E.V*C.oL)@f/E.mass))
 d.update(steps=steps,rejections=rej,t=t,dt=dt,history=history,repair=report(f,out),raw=report(fr,outr),distribution_L1=float(E.V@abs(f-fr)/E.mass),cumulative_repair_L1=sum(x['l1'] for x in history),cumulative_raw_error=sum(x['error'] for x in history),max_relative_repair=max(x['relative_l1'] for x in history),max_alpha_over_dt=max(x['alpha'] for x in history),max_capacity_frac=max(x['max_frac'] for x in history),max_step_energy_residual=max(abs(x['after']) for x in history));return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=1);p.add_argument('--nphase',type=int,default=32);p.add_argument('--nrad',type=int,default=10);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=4);p.add_argument('--rout',type=float,default=2.);p.add_argument('--nz',type=int,default=64);p.add_argument('--t',type=float,default=.1);p.add_argument('--dt',type=float,default=.005);p.add_argument('--tilt',action='store_true');a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt)))
