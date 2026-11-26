"""Frozen reference orbits with evolving nonlinear collisions + linearized
Vlasov-Poisson feedback on the initial lowered Bose equilibrium.
This omits delta-force acting on delta-f and does not claim nonlinear gravity.
"""
import argparse,json
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix
from numpy.polynomial.legendre import leggauss
from lowered_bose_halo_orbit import build,Halo,TIME,open_shift,shift
from lowered_bose_orbit_tilt import Tilt
from lowered_bose_orbit_batch import Streamer

class Feedback:
 def __init__(self,H,C,R,ng=8):
    self.C=C;self.H=H;self.R=R;self.ns=len(R)-1;self.nphase=C.nphase
    j=C.shell;J=coo_matrix((np.ones(len(j)),(j,np.arange(len(j)))),shape=(self.ns,len(j))).tocsr();self.W=(J@C.W).tocsr();self.initial=self.W@C.f0
    rows=[];cols=[];vals=[];x,w=leggauss(ng)
    for O,ed,base in C.orbits:
        if O.s<=0:continue
        d=np.diff(ed);s=O.s;f=1/np.expm1(H.eta-s);fp=f*(1+f)-(H.d if H.model=='quadratic' else 0.)
        for j,(a,b) in enumerate(zip(R[:-1],R[1:])):
            for l,h in O.intervals(a,b):
                ii=np.arange(max(0,np.searchsorted(ed,l,side='right')-1),min(len(d),np.searchsorted(ed,h,side='left')));aa=np.maximum(l,ed[ii]);bb=np.minimum(h,ed[ii+1]);keep=bb>aa;ii=ii[keep];aa=aa[keep];bb=bb[keep]
                if not len(ii):continue
                th=(aa[:,None]+bb[:,None])/2+(bb-aa)[:,None]*x/2;r,p=O.coords(th);frac=np.clip((r**3-a**3)/(b**3-a**3),0,1);wt=-fp*np.sqrt(2)*p/r**2*(bb-aa)[:,None]*w/2/d[ii,None]/TIME
                left=(wt*(1-frac)).sum(1);right=(wt*frac).sum(1)
                rows.extend(base+ii);cols.extend([j]*len(ii));vals.extend(left)
                rows.extend(base+ii);cols.extend([j+1]*len(ii));vals.extend(right)
    G=coo_matrix((vals,(rows,cols)),shape=(len(C.V),self.ns+1)).tocsr()
    rev=np.arange(len(C.V)).reshape(-1,self.nphase)[:,::-1].ravel();self.G=(G-G[rev])/2
    self.antisym_raw=float(abs(G+G[rev]).max());self.number_cols=float(abs(self.G.T@C.V).max());self.energy_cols=float(abs(self.G.T@(C.V*C.e)).max());self.q0=H.sol.sol(np.clip(R,H.eps,H.rt))[1];self.q0=np.where(R<H.eps,R**3/3,self.q0);self.q0=np.where(R>H.rt,H.Qt,self.q0)
 def enclosed(self,f):return np.r_[0,np.cumsum(self.W@f-self.initial)/self.H.Ic]
 def deriv(self,f):return self.G@self.enclosed(f)
 def report(self,f):
    q=self.enclosed(f);r=self.R;dg=self.deriv(f);rel=np.divide(q,self.q0,out=np.zeros_like(q),where=self.q0>0);dM=self.W@f-self.initial
    x,w=leggauss(12);a=r[:-1];b=r[1:];rr=(a[:,None]+b[:,None])/2+(b-a)[:,None]*x/2;fr=(rr**3-a[:,None]**3)/(b[:,None]**3-a[:,None]**3);qq=q[:-1,None]+fr*np.diff(q)[:,None];field=-self.H.Ic/2*(np.sum((b-a)[:,None]*w/2*qq*qq/rr**2)+q[-1]**2/r[-1]);mass=self.C.V@self.C.f0
    return dict(q_edges=q.tolist(),max_abs_force=float(np.max(abs(q[1:]/r[1:]**2))),max_relative_force=float(abs(rel[1:]).max()),relative_outside_point1=float(abs(rel[r>=.1]).max()),max_source=float(abs(dg).max()),force_Ndot=float(self.C.V@dg),force_Edot=float((self.C.V*self.C.e)@dg),perturbation_field_energy=float(field/mass),shell_mass_residual=float(dM.sum()-(self.C.V@(f-self.C.f0))))

def run(H,M=1,K=4,nL=3,nphase=64,nrad=12,m=5,kphase=6,rout=2.,nz=128,tmax=.1,dt=.005,tilt=False):
    C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz);P=Tilt(C) if tilt else C;V=C.V;e=C.e;mass=V@C.f0;F=Feedback(H,C,C.R);f0=C.f0.copy();f=f0.copy();fr=f0.copy();out=np.zeros(2);outr=np.zeros(2);t=0.;steps=0;rej=0
    d.update(dt=dt,tmax=tmax,tilt=tilt,feedback_matrix_nonzero=F.G.nnz,feedback_antisym_raw=F.antisym_raw,force_N_col=F.number_cols,force_E_col=F.energy_cols,radial_edges=C.R.tolist())
    S=Streamer(C,TIME)
    def stream(f,h):
        z,out,_=S.apply(f,h);return z,out[:2]
    def step(f,h,gravity):
        x,o=stream(f,h/2);g=P.derivative(x)+(F.deriv(x) if gravity else 0.);p=x+h*g
        if p.min()<0:return None,None
        g2=P.derivative(p)+(F.deriv(p) if gravity else 0.);y=.5*(x+p+h*g2)
        if y.min()<0:return None,None
        y,q=stream(y,h/2);return y,o+q
    def rep(f,t,out):
        frow=C.stats(f);frow.update(t=t,gaussian=(C.dens.T@(V*f)/mass).tolist(),gaussian_change=(C.dens.T@(V*(f-f0))/mass).tolist(),cutoff=float(V[e>0]@f[e>0]/mass),unbound=float(V[e>C.e0]@f[e>C.e0]/mass),outside=float((V*C.outside)@f/mass),count_ledger=float((V@f+out[0])/mass-1),reference_energy_ledger=float(((V*e)@f+out[1]-(V*e)@f0)/mass),out_count=float(out[0]/mass),**F.report(f));return frow
    rows=[dict(vp=rep(f,0,out),frozen=rep(fr,0,outr))]
    while t<tmax-1e-12:
        h=min(dt,tmax-t)
        while True:
            a,b=step(f,h,True);c,dd=step(fr,h,False)
            if a is not None and c is not None:break
            h*=.5;rej+=1
            if h<1e-9:raise RuntimeError('positive VP step too small; linearization or time step failure')
        f,fr=a,c;out+=b;outr+=dd;t+=h;steps+=1
    rows.append(dict(vp=rep(f,t,out),frozen=rep(fr,t,outr)));d.update(steps=steps,rejections=rej,rows=rows);return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=3);p.add_argument('--nphase',type=int,default=64);p.add_argument('--nrad',type=int,default=12);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=6);p.add_argument('--rout',type=float,default=2);p.add_argument('--nz',type=int,default=128);p.add_argument('--t',type=float,default=.1);p.add_argument('--dt',type=float,default=.005);p.add_argument('--tilt',action='store_true');a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt)))
