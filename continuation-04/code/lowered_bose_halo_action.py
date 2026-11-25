"""Nonlinear spherical residual Hamiltonian transport on common L² orbit grid.
Full force on f, reference-orbit streaming, nonlinear C22. Finite-domain prototype.
"""
import argparse,json
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from scipy.sparse import coo_matrix
from lowered_bose_halo_orbit import build,Halo,Orbit,TIME
from lowered_bose_halo_vp import Feedback
from lowered_bose_orbit_batch import Streamer
from lowered_bose_orbit_tilt import Tilt

def mm(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)

class Field:
 def __init__(self,H,C,nz=64):
    assert C.common
    self.H=H;self.C=C;self.n=C.nphase;self.no=len(C.orbits);self.r=C.R;self.ns=len(self.r)-1;self.e=np.array([-o.s for o,_,_ in C.orbits]);self.L=np.array([o.L for o,_,_ in C.orbits]);self.open=np.array([o.open for o,_,_ in C.orbits]);self.w=C.Lweights/(8*np.sqrt(2)*TIME);self.edges=C.orbits[0][1];assert all(np.array_equal(ed,self.edges) for _,ed,_ in C.orbits)
    self.d=np.diff(self.edges);self.tc=(self.edges[:-1]+self.edges[1:])/2;self.tl=np.r_[self.tc[-1]-1,self.tc[:-1]];self.tr=np.r_[self.tc[1:],self.tc[0]+1]
    j=C.shell;J=coo_matrix((np.ones(len(j)),(j,np.arange(len(j)))),shape=(self.ns,len(j))).tocsr();self.W=(J@C.W).tocsr();self.initial=self.W@C.f0
    self.lo=np.zeros(self.no,int);self.hi=np.zeros(self.no,int);self.prev=np.arange(self.no);self.next=np.arange(self.no);self.el=np.zeros(self.no);self.er=np.zeros(self.no);rad=[];lf=[];rf=[];wf=[];ef=[];Lf=[];tops=[];psi=C.psi;R=C.rmax
    order=np.lexsort((self.e,self.L));cuts=np.r_[0,np.flatnonzero(np.diff(self.L[order]))+1,len(order)]
    for a,b in zip(cuts[:-1],cuts[1:]):
        idx=order[a:b];L=self.L[idx[0]];es=self.e[idx];assert np.allclose(np.diff(es),C.de,rtol=1e-10,atol=1e-12);assert np.max(abs(self.w[idx]-self.w[idx[0]]))<1e-14
        emin=es[0]+C.de/2-C.energy_width[idx[0]];eface=np.r_[emin,(es[:-1]+es[1:])/2,es[-1]+C.de/2]
        k=len(rad);self.lo[idx]=k+np.arange(len(idx));self.hi[idx]=self.lo[idx]+1;self.prev[idx[1:]]=idx[:-1];self.next[idx[:-1]]=idx[1:];self.el[idx]=es-eface[:-1];self.er[idx]=eface[1:]-es
        fit=minimize_scalar(lambda r:-float(psi(r))+L/(2*r*r),bounds=(R*1e-12,R),method='bounded',options={'xatol':1e-10});rc0=R if -float(psi(R))+L/(2*R*R)<fit.fun else fit.x
        rad.append(np.full(self.n+1,rc0));lf.append(-1);rf.append(idx[0]);wf.append(self.w[idx[0]]);ef.append(eface[0]);Lf.append(L)
        for j,ee in enumerate(eface[1:]):
            s=-ee;rs=brentq(lambda r:float(psi(r))-s,0,R,xtol=1e-12) if s>float(psi(R)) else R
            g=lambda r:2*r*r*(float(psi(r))-s);sol=minimize_scalar(lambda r:-g(r),bounds=(0,rs),method='bounded',options={'xatol':1e-9});rc=sol.x;gm=-sol.fun
            if g(rs)>gm:rc=rs;gm=g(rs)
            if gm<=L:raise RuntimeError(('inaccessible energy face',L,ee,gm))
            isopen=bool(rs==R and g(R)>L);o=Orbit(psi,s,L,rs,rc,isopen,nz,middle=True);r,p=o.coords(self.edges);r[-1]=r[0];rad.append(r);lf.append(idx[j]);rf.append(idx[j+1] if j+1<len(idx) else -1);wf.append(self.w[idx[j]]);ef.append(ee);Lf.append(L)
        tops.append(len(rad)-1)
    self.rad=np.array(rad);self.lf=np.array(lf);self.rf=np.array(rf);self.wf=np.array(wf);self.ef=np.array(ef);self.Lf=np.array(Lf);self.tops=np.array(tops);self.ngroups=len(tops);self.valid_l=self.lf>=0;self.valid_r=self.rf>=0
    r=np.clip(self.rad,0,C.rmax);j=np.clip(np.searchsorted(self.r,r,side='right')-1,0,self.ns-1);a=self.r[j];b=self.r[j+1];self.vj=j;self.vI=np.divide(1,r,out=np.zeros_like(r),where=r>0)-1/b;self.vJ=(b*b-r*r)/2
    self.shell_a=self.r[:-1];self.shell_b=self.r[1:];self.shell_vol=self.shell_b**3-self.shell_a**3
    self.q0=H.sol.sol(np.clip(self.r,H.eps,H.rt))[1];self.q0=np.where(self.r<H.eps,self.r**3/3,self.q0);self.q0=np.where(self.r>H.rt,H.Qt,self.q0)
    self.F=Feedback(H,C,C.R)
 def enclosed(self,f):return np.r_[0,np.cumsum((self.W@f-self.initial)/self.H.Ic)]
 def potential(self,q):
    a,b=self.shell_a,self.shell_b;c1=np.diff(q)/self.shell_vol;c0=q[:-1]-c1*a**3;ia=np.divide(1,a,out=np.zeros_like(a),where=a>0)-1/b;I=c0*ia+c1*(b*b-a*a)/2;suf=np.r_[np.cumsum(I[::-1])[::-1],0.];j=self.vj
    return -q[-1]/self.r[-1]-suf[j+1]-c0[j]*self.vI-c1[j]*self.vJ
 def fluxes(self,q):
    h=self.potential(q);h[:,-1]=h[:,0];Fe=-np.diff(h,axis=1)*self.wf[:,None];Ft=(h[self.hi]-h[self.lo])*self.w[:,None]
    div=Fe[self.hi]-Fe[self.lo]+np.diff(Ft,axis=1);return Fe,Ft,float(abs(div).max())
 def faces(self,f,q=None,order=2,signed=False):
    if q is None:q=self.enclosed(f)
    Fe,Ft,div=self.fluxes(q);z=f.reshape(self.no,self.n);es=np.zeros_like(z);ts=np.zeros_like(z)
    if order==2:
        dl=self.e-self.e[self.prev];dr=self.e[self.next]-self.e;l=np.divide(z-z[self.prev],dl[:,None],out=np.zeros_like(z),where=dl[:,None]>0);r=np.divide(z[self.next]-z,dr[:,None],out=np.zeros_like(z),where=dr[:,None]>0);es=mm(l,r)
        ts=mm((z-np.roll(z,1,axis=1))/(self.tc-self.tl),(np.roll(z,-1,axis=1)-z)/(self.tr-self.tc));ts[self.open,0]=ts[self.open,-1]=0.
        if not signed:
            es=np.clip(es,-z/self.er[:,None],z/self.el[:,None]);ts=np.clip(ts,-2*z/self.d,2*z/self.d)
    ezl=z+es*self.er[:,None];ezr=z-es*self.el[:,None];a=np.zeros_like(Fe);b=np.zeros_like(Fe);a[self.valid_l]=ezl[self.lf[self.valid_l]];b[self.valid_r]=ezr[self.rf[self.valid_r]];EF=Fe*np.where(Fe>=0,a,b)
    tl=np.hstack((z[:,-1:]+ts[:,-1:]*self.d[-1]/2,z+ts*self.d/2));tr=np.hstack((z-ts*self.d/2,z[:,0:1]-ts[:,0:1]*self.d[0]/2));tl[self.open,0]=0.;tr[self.open,-1]=0.;TF=Ft*np.where(Ft>=0,tl,tr)
    return EF,TF,Fe,Ft,div
 def counts(self,EF,TF):
    g=-(EF[self.hi]-EF[self.lo]+np.diff(TF,axis=1));oe=EF[self.tops].sum(1);ot=(TF[self.open,-1]-TF[self.open,0]);out=np.array([oe.sum()+ot.sum(),oe@self.ef[self.tops]+ot@self.e[self.open],oe@self.Lf[self.tops]+ot@self.L[self.open]])
    return g,out,oe,ot
 def rhs(self,f,q=None,order=2,diagnostic=False,signed=False):
    EF,TF,Fe,Ft,div=self.faces(f,q,order,signed);g,out,oe,ot=self.counts(EF,TF)
    if diagnostic:
        outgoing=np.maximum(Fe[self.hi],0)-np.minimum(Fe[self.lo],0)+np.maximum(Ft[:,1:],0)-np.minimum(Ft[:,:-1],0);cfl=.35/max(float((outgoing.ravel()/self.C.V).max()),1e-300)
        return g.ravel()/self.C.V,out,dict(flux_div=div,cfl=cfl,Ndot=float(g.sum()+out[0]),Ldot=float(self.L@g.sum(1)+out[2]),Edot=float(self.e@g.sum(1)+out[1]),energy_boundary=float(oe.sum()),spatial_boundary=float(ot.sum()))
    return g.ravel()/self.C.V,out
 def field_energy(self,f):return self.F.report(f)['perturbation_field_energy']
 def control(self):
    q=.001*self.q0;f=self.C.f0;g,o,d=self.rhs(f,q,2,True);lin=self.F.G@q;V=self.C.V;rev=g.reshape(self.no,self.n)[:,::-1].ravel();odd=(g-rev)/2;even=(g+rev)/2;norm=V@abs(lin)
    return dict(**d,linear_L1=float(norm),full_relative_L1=float(V@abs(g-lin)/norm),odd_relative_L1=float(V@abs(odd-lin)/norm),even_relative_L1=float(V@abs(even)/norm),density_source_L1=float(abs(self.W@g).sum()/norm),linear_density_source_L1=float(abs(self.W@lin).sum()/norm),exact_zero=float(abs(self.rhs(f,np.zeros_like(q))[0]).max()))

def run(H,M=1,K=4,nL=1,nphase=32,nrad=10,m=5,kphase=4,rout=2.,nz=64,tmax=.1,dt=.005,tilt=False,deltaf=False,order=2):
    C,d=build(H,M,K,nL,nphase,nrad,m,kphase,rout,nz,True);P=Tilt(C) if tilt else C;F=Field(H,C,nz);S=Streamer(C,TIME);V=C.V;e=C.e;f0=C.f0.copy();f=f0.copy();fr=f0.copy();mass=V@f0;out=np.zeros(2);outr=np.zeros(2);steps=rej=0;t=0.;g0=C.derivative(f0);control=F.control();d.update(control=control,action_L_groups=F.ngroups,action_faces=len(F.rad),initial_source=[float(V[e>0]@g0[e>0]/mass),float(V[e>C.e0]@g0[e>C.e0]/mass)],deltaf=deltaf,order=order,dt=dt,tmax=tmax)
    def gravity(f):
        if not deltaf:return F.rhs(f,order=order)
        q=F.enclosed(f);g,o=F.rhs(f-f0,q,order,signed=True);return g+F.F.G@q,o
    def step(f,h,gr):
        x,o,_=S.apply(f,h/2);a,oa=gravity(x) if gr else (0.,np.zeros(3));g=P.derivative(x)+a;p=x+h*g
        if p.min()<-1e-13:return None,None
        p=np.maximum(p,0);a2,ob=gravity(p) if gr else (0.,np.zeros(3));y=.5*(x+p+h*(P.derivative(p)+a2))
        if y.min()<-1e-13:return None,None
        y=np.maximum(y,0);z,oo,_=S.apply(y,h/2);return z,o[:2]+oo[:2]+h*(oa[:2]+ob[:2])/2
    def report(f,t,out):
        field=F.F.report(f);E=float(((V*e)@(f-f0)+out[1])/mass);g,o,gd=F.rhs(f,order=order,diagnostic=True);return dict(t=t,Nledger=float((V@f+out[0])/mass-1),reference_energy_change=E,field_energy=field['perturbation_field_energy'],total_energy_error=E+field['perturbation_field_energy'],min=float(f.min()),max_force_rel=field['max_relative_force'],gaussian_change=(C.dens.T@(V*(f-f0))/mass).tolist(),cutoff=float(V[e>0]@f[e>0]/mass),unbound=float(V[e>C.e0]@f[e>C.e0]/mass),outside=float((V*C.outside)@f/mass),out_count=float(out[0]/mass),gravity_ledger=gd)
    rows=[dict(full=report(f,0,out),frozen=report(fr,0,outr))]
    while t<tmax-1e-12:
        h=min(dt,tmax-t)
        while True:
            a,b=step(f,h,True);c,dd=step(fr,h,False)
            if a is not None and c is not None:break
            h*=.5;rej+=1
            if h<1e-8:raise RuntimeError('positivity failure, dt below floor')
        f,fr=a,c;out+=b;outr+=dd;t+=h;steps+=1
    rows.append(dict(full=report(f,t,out),frozen=report(fr,t,outr)));d.update(steps=steps,rejections=rej,rows=rows);return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--model',default='linear');p.add_argument('--delta',type=float,default=.1);p.add_argument('--M',type=int,default=1);p.add_argument('--K',type=int,default=4);p.add_argument('--nL',type=int,default=1);p.add_argument('--nphase',type=int,default=32);p.add_argument('--nrad',type=int,default=10);p.add_argument('--m',type=int,default=5);p.add_argument('--kphase',type=float,default=4);p.add_argument('--rout',type=float,default=2);p.add_argument('--nz',type=int,default=64);p.add_argument('--t',type=float,default=.1);p.add_argument('--dt',type=float,default=.005);p.add_argument('--tilt',action='store_true');p.add_argument('--deltaf',action='store_true');p.add_argument('--order',type=int,default=2);a=p.parse_args();print(json.dumps(run(Halo(delta=a.delta,model=a.model),a.M,a.K,a.nL,a.nphase,a.nrad,a.m,a.kphase,a.rout,a.nz,a.t,a.dt,a.tilt,a.deltaf,a.order)))
