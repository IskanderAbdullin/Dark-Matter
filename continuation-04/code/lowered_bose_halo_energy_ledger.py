"""Exact quadratic discrete shell-field functional and split work ledger."""
import argparse,json,numpy as np
from numpy.polynomial.legendre import leggauss
from lowered_bose_halo_action_limiter import LimitedField,build,Halo,TIME,Streamer,Tilt

class Energy:
 def __init__(self,F,n=12):
    self.F=F;self.V=F.C.V;self.e=F.C.e;self.mass=self.V@F.C.f0;self.W=F.W;self.initial=F.initial
    a,b=F.r[:-1],F.r[1:];x,w=leggauss(n);rr=(a[:,None]+b[:,None])/2+(b-a)[:,None]*x/2;ww=(b-a)[:,None]*w/2/rr**2
    ns=len(a);D=np.broadcast_to(np.arange(ns)[None,None,:]<np.arange(ns)[:,None,None],(ns,n,ns)).astype(float).copy();D[np.arange(ns),:,np.arange(ns)]=(rr**3-a[:,None]**3)/(b[:,None]**3-a[:,None]**3);D=D.reshape(-1,ns);self.A=(D.T@(ww.ravel()[:,None]*D)+1/F.r[-1])/F.H.Ic
 def shell(self,f):return self.W@f-self.initial
 def field(self,f):
    m=self.shell(f);return -.5*m@self.A@m
 def hamiltonian(self,f):return self.e-(self.W.T@(self.A@self.shell(f)))/self.V
 def difference(self,f,g,out=0.):
    a,b=self.shell(f),self.shell(g);phi=-.5*self.A@(a+b);return np.array([(self.V*self.e)@(g-f)+out,phi@(b-a)])/self.mass
 def split(self,f,g,dc,og):
    a,b=self.shell(f),self.shell(g);phi=-.5*self.A@(a+b);dg=g-f-dc
    return np.array([[(self.V*self.e)@dc,phi@(self.W@dc)],[(self.V*self.e)@dg+og,phi@(self.W@dg)]])/self.mass

def run(H,M=1,K=4,nL=1,nphase=32,nrad=10,m=5,kphase=4,rout=2.,nz=64,tmax=.1,dt=.005,tilt=False,mode='limiter'):
 C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz,True);P=Tilt(C) if tilt else C;F=LimitedField(H,C,nz) if mode=='limiter' else None;
 if F is None:
    from lowered_bose_halo_action_variants import WorkField
    F=WorkField(H,C,nz,mode)
 S=Streamer(C,TIME);E=Energy(F);f=C.f0.copy();out=np.zeros(2);t=0.;steps=rej=0;ledger=np.zeros((3,2));absledger=np.zeros((3,2));maxcheck=0.;max_shell=0.;minstep=dt
 while t<tmax-1e-12:
    h=min(dt,tmax-t)
    while True:
        x,o,_=S.apply(f,h/2);c=P.derivative(x);p,oa,da=F.euler(x,c,h)
        if p is not None:
            c2=P.derivative(p);y,ob,db=F.euler(p,c2,h)
            if y is not None:break
        rej+=1;h*=.5
        if h<1e-8:raise RuntimeError('positivity failure')
    y=(x+y)/2;og=(oa+ob)/2;z,oo,_=S.apply(y,h/2);dc=h*(c+c2)/2
    q=np.vstack((E.difference(f,x,o[1])+E.difference(y,z,oo[1]),E.split(x,y,dc,og[1])))
    maxcheck=max(maxcheck,float(abs(q.sum(0)-E.difference(f,z,o[1]+oo[1]+og[1])).max()));max_shell=max(max_shell,float(abs(E.W@dc).sum()/E.mass));ledger+=q;absledger+=abs(q);out+=o[:2]+oo[:2]+og[:2];f=z;t+=h;steps+=1;minstep=min(minstep,h)
 field=E.field(f)/E.mass;field_ref=F.F.report(f)['perturbation_field_energy'];direct=E.difference(C.f0,f,out[1]);ham=E.hamiltonian(f);rng=np.random.default_rng(42);g=rng.normal(size=len(f))*f;eps=1e-5;fd=(E.field(f+eps*g)-E.field(f-eps*g))/(2*eps);ana=(E.V*(ham-E.e))@g
 d.update(steps=steps,rejections=rej,minstep=minstep,field=field,field_reference=field_ref,gradient_fd=fd,gradient_exact=float(ana),gradient_error=float(fd-ana),ledger=ledger.tolist(),abs_ledger=absledger.tolist(),labels=['reference stream','C22','residual gravity'],direct=direct.tolist(),ledger_residual=(ledger.sum(0)-direct).tolist(),max_step_closure=maxcheck,max_collision_shell_mass_change=max_shell,Nledger=float((E.V@f+out[0])/E.mass-1),out_count=float(out[0]/E.mass),gaussian_change=(C.dens.T@(E.V*(f-C.f0))/E.mass).tolist());return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=1);p.add_argument('--nphase',type=int,default=32);p.add_argument('--nrad',type=int,default=10);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=4);p.add_argument('--rout',type=float,default=2.);p.add_argument('--nz',type=int,default=64);p.add_argument('--t',type=float,default=.1);p.add_argument('--dt',type=float,default=.005);p.add_argument('--tilt',action='store_true');p.add_argument('--mode',choices=['limiter','parity','work','adjoint'],default='limiter');a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt,a.mode)))
