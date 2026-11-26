"""Short closed-support joint resonant evolution with number/kinetic-energy cells."""
import argparse,json
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.optimize import brentq
from joint_nonlinear_prompt_return_moments_reference import geometry
from fermion_two_moment_closure_reference import infer as fd
from regular_bose_two_moment_closure_reference import infer as be
from sm_dark_entropy_exchange_reference import StandardModel,MPL
from three_species_decay_onsager_reference import M,m,mp

@lru_cache(None)
def rule(n):
    return roots_legendre(n)
def kin(p,mass):
    return p*p/(np.hypot(p,mass)+mass)

class Cells:
    def __init__(self,p,mass,g):
        self.p=np.asarray(p);self.mass=mass;self.g=g
        self.low=kin(p[:-1],mass);self.hi=kin(p[1:],mass);self.dk=self.hi-self.low
        self.vol=g*np.diff(p**3)/(6*np.pi**2)
        a,b=rule(96)
        q=(p[:-1,None]*(1-a)+p[1:,None]*(1+a))/2
        w=q*q*b
        self.w=w/w.sum(axis=1)[:,None]
        self.u=(kin(q,mass)-self.low[:,None])/self.dk[:,None]
        self.mean=self.low+self.dk*np.sum(self.w*self.u,axis=1)

    def f(self,ab,e,idx):
        a,b=ab
        u=np.clip(((e-self.mass)-self.low[idx])/self.dk[idx],0.,1.)
        s=a[idx]+b[idx]*u
        if self.g==4:return expit(-s)
        a=np.exp(-s)
        return a/(-np.expm1(-s))

    def reconstruct(self,state):
        n,k=state
        if np.any(n<0.) or np.any((n<=0.)&(k!=0.)): raise ValueError("count")
        fill=n/self.vol
        if self.g==4 and np.any(fill>=1.):raise ValueError("Pauli")
        ar,br=np.full(len(n),1000.),np.zeros(len(n))
        for i in np.flatnonzero(n>0.):
            x=(k[i]/n[i]-self.low[i])/self.dk[i]
            if not 0.<x<1.:raise ValueError(("support",i,x))
            ref=float(self.w[i]@self.u[i])
            if abs(ref-x)<5e-14:
                A,B=(np.log1p(fill[i])-np.log(fill[i]),0.) if self.g==1 else (np.log1p(-fill[i])-np.log(fill[i]),0.)
            elif self.g==1:
                A,B=be(*self.p[i:i+2],self.mass,fill[i],x)
            else:
                A,B=fd(*self.p[i:i+2],self.mass,fill[i],x)
                A=-A
            ar[i],br[i]=A,B
        s=ar[:,None]+br[:,None]*self.u
        if self.g==1:
            if np.min(s)<=0.:raise ValueError("BE domain")
            f=np.exp(-s)/(-np.expm1(-s))
        else:
            f=expit(-s)
        np1=np.sum(self.w*f,axis=1)*self.vol
        kp=np.sum(self.w*f*(self.low[:,None]+self.dk[:,None]*self.u),axis=1)*self.vol
        take=n>1e-240
        err=max(np.max(abs(np1[take]/n[take]-1)),np.max(abs(kp[take]/k[take]-1)))
        if err>3e-7:raise ValueError(("moment recovery",err))
        entropy=s*f
        if self.g==1:
            entropy-=np.log(-np.expm1(-s))
        else:
            entropy+=np.logaddexp(0.,-s)
        return (ar,br),float(err),float(np.sum(self.vol[:,None]*self.w*entropy))

class Collision:
    def __init__(self,p,order=5,angle=4):
        self.c=[Cells(p,mp,1),Cells(p,m,4)]
        parent,w,mid,x,ip,iq=geometry(p,p,2*p[-1]+M,order)
        a,v=rule(angle)
        active=(x>0)&(ip<len(p)-1)&(iq<len(p)-1)
        if p[0]>0.:
            # Geometry includes the full parent decay interval; clipped
            # integer indices alone do not enforce a nonzero inner cut.
            active&=(mid>=np.hypot(p[0],mp))&((parent[:,None]-mid)>=np.hypot(p[0],m))
        y,j=np.nonzero(active)
        self.E=parent
        self.w=w
        self.y=np.repeat(y,len(a))
        self.n=len(parent)
        self.I=np.repeat(ip[y,j],len(a))
        self.J=np.repeat(iq[y,j],len(a))
        self.ph=(mid[y,j,None]+x[y,j,None]*a/2).ravel()
        self.ch=self.E[self.y]-self.ph
        self.pre=(x[y,j,None]*v/2).ravel()
        self.fac=w[self.y]
        self.k1=self.ph-mp
        self.k2=self.ch-m

    def rate(self,shapes):
        f=self.c[0].f(shapes[0],self.ph,self.I)
        g=self.c[1].f(shapes[1],self.ch,self.J)
        A=self.pre*f*g
        B=self.pre*(1+f)*(1-g)
        x=np.bincount(self.y,weights=A,minlength=self.n)
        y=np.bincount(self.y,weights=B,minlength=self.n)
        z=x+y
        # Avoid subtracting almost identical O(1) parent blocking factors.
        zz=np.divide(self.fac,z[self.y],out=np.zeros_like(self.fac),where=z[self.y]>0.)
        F=zz*(B*x[self.y]-A*y[self.y])
        n=len(self.c[0].vol)
        D=np.array([[np.bincount(k,weights=F,minlength=n),np.bincount(k,weights=F*E,minlength=n)]
                         for k,E in ((self.I,self.k1),(self.J,self.k2))])
        return D

    def entropy_control(self,T):
        shapes=[]
        for c in self.c:
            mu=-.23*T if c.g==1 else -5.7*T
            shapes.append(np.array([(c.mass+c.low-mu)/T,c.dk/T]))
        eq=self.rate(shapes)
        f=self.c[0].f(shapes[0],self.ph,self.I)
        g=self.c[1].f(shapes[1],self.ch,self.J)
        gross0=float(np.sum(self.fac*self.pre*f*g))
        equilibrium=np.max(abs(eq[:,0,:]))/gross0
        assert equilibrium<4e-13
        for j in range(2):
            shapes[j][0]+=0.09*np.sin(np.arange(shapes[j].shape[1])*.8+j)
            shapes[j][1]*=(.92+.36*j)
        net=self.rate(shapes)
        terms=[]
        for i,c in enumerate(self.c):
            a,b=shapes[i]
            terms.append(float((a-b*c.low/c.dk)@net[i,0]+(b/c.dk)@net[i,1]))
        S=sum(terms)
        f=self.c[0].f(shapes[0],self.ph,self.I)
        g=self.c[1].f(shapes[1],self.ch,self.J)
        a=self.pre*f*g
        b=self.pre*(1+f)*(1-g)
        u=np.bincount(self.y,weights=a,minlength=self.n)[self.y]
        v=np.bincount(self.y,weights=b,minlength=self.n)[self.y]
        mask=(a>0)&(b>0)&(u>0)&(v>0)
        A=a[mask]*v[mask];B=b[mask]*u[mask]
        S2=np.sum(self.fac[mask]*(B-A)/(u[mask]+v[mask])*(np.log(B)-np.log(A)))
        error=abs(S2/S-1)
        assert S>=0 and S2>=0 and error<2e-11,(S,S2,error)
        return dict(exact_continuum_KMS_from_cell_shapes=equilibrium,
                    microscopic_entropy_relative_check=error,entropy_production=S)


def maxent(c,N,Et,guess):
    # Separate numbers; common temperature chosen by true total energy.
    # All energies and populations are integrated on the same finite box.
    def evaluate(T,cs,n):
        E=cs.mass+cs.low[:,None]+cs.dk[:,None]*cs.u
        weight=cs.w*cs.vol[:,None]
        def occupation(eta):
            z=E/T-eta
            return expit(-z) if cs.g==4 else np.exp(-z)/(-np.expm1(-z))
        hi=50. if cs.g==4 else cs.mass/T-1e-10
        eta=brentq(lambda z:np.sum(weight*occupation(z))/n-1.,-150.,hi,xtol=2e-13)
        f=occupation(eta);z=E/T-eta
        S=z*f+np.logaddexp(0.,-z) if cs.g==4 else z*f-np.log(-np.expm1(-z))
        return (float(np.sum(weight*f*E)),float(np.sum(weight*S)),np.sum(weight*f,axis=1),eta)
    def state(T):return [evaluate(T,j,n) for j,n in zip(c,N)]
    tc=brentq(lambda T:sum(j[0] for j in state(T))/Et-1.,.25*guess,4*guess,xtol=1e-14)
    return tc,state(tc)


def run(row,eos,step,final,dt,method="heun",order=5,verify=False,angle=4,cap=0.):
    T=row["TF"];s=eos.state(T)[2];H=np.sqrt(eos.state(T)[0]/3)/MPL
    pp=np.array(row["scalar_differential_bins"])*T
    last=len(pp)-1 if cap<=0. else np.searchsorted(pp,T*cap,side="right")-1
    ids=np.unique(np.r_[np.arange(0,last+1,step),last])
    ker=Collision(pp[ids],order,angle)
    sp=np.array(row["scalar_normalized_counts"])
    fp=np.array(row["stable_normalized_counts"])
    n=np.array([s*row["stable_yields"][1]*np.add.reduceat(sp[:last],ids[:-1]),
                s*sum(row["stable_yields"])*np.add.reduceat(fp[:last],ids[:-1])])
    y=np.array([n,n*np.array([c.mean for c in ker.c])]).transpose(1,0,2)
    norm=y.sum()
    U=y.copy()/norm
    def reconstruction(a):
        r=[c.reconstruct(z*norm) for c,z in zip(ker.c,a)]
        return [q[0] for q in r],max(q[1] for q in r),sum(q[2] for q in r)
    def rhs(a,shape):return ker.rate(shape)/(H*norm)
    v,err,S0=reconstruction(U)
    x=0.;nstep=reject=0;du=0;S=S0
    controls=ker.entropy_control(T) if verify else None
    conserved=U[:,0,:].sum(axis=1)
    def energy(a):return a[:,1,:].sum()+a[0,0,:].sum()*mp+a[1,0,:].sum()*m
    E0=energy(U)
    h_history=[]
    while x<final:
        h=min(dt,final-x)
        D=rhs(U,v)
        exact=None
        for jj in range(14):
            try:
                C=U+h*D
                z,e,SS=reconstruction(C)
                if method=="heun":
                    C=U+.5*h*(D+rhs(C,z))
                    z,e,SS=reconstruction(C)
                if SS<S-3e-8*S0:
                    raise ValueError("entropy decreased")
                exact=(C,z,e,SS)
                break
            except ValueError:
                h*=.25;reject+=1
        if exact is None:raise RuntimeError("no conservative admissible step")
        U,v,e,S=exact
        x+=h;err=max(err,e);nstep+=1
        h_history.append(h)
        if nstep>1300:raise RuntimeError("explicit time boundary")
    inv=max(abs(U[:,0,:].sum(axis=1)/conserved-1))
    energy_error=abs(energy(U)/E0-1)
    assert max(inv,energy_error)<2e-11
    qe,es=maxent(ker.c,conserved*norm,E0*norm,T)
    Smax=sum(j[1] for j in es)
    assert Smax>S0 and S<=Smax*(1+2e-8)
    ideal_scalar_heat=(es[0][0]-(y[0,1,:].sum()+mp*y[0,0,:].sum()))/(E0*norm)
    L1=[float(np.sum(abs(U[i,0,:]*norm-es[i][2]))/(conserved[i]*norm)) for i in range(2)]
    return dict(T=T,bins=len(ids)-1,method=method,log_time_H_units=final,entropy_controls=controls,
       true_finite_box_common_temperature=qe,
       physical_cutoff_p_over_T=float(pp[last]/T),
       removed_by_initial_support_fraction=[float(sp[last:].sum()),float(fp[last:].sum())],
       parent_and_angle_orders=[order,angle],
       finite_box_common_temperature_limit_scalar_heat=ideal_scalar_heat,
       fractional_progress_towards_maximum_entropy=(S-S0)/(Smax-S0),
       terminal_total_variation_from_common_temperature=[j/2 for j in L1],
       steps=nstep,rejections=reject,minimum_accepted_step_H_units=min(h_history),
       total_energy_relative_error=energy_error,separate_number_relative_error=inv,
       maximum_moment_shape_error=err,relative_entropy_gain=float(S/S0-1),
       fractional_scalar_heat_of_initial_dark_energy=float((U[0,1,:].sum()+mp*U[0,0,:].sum()-energy(np.r_[y[:1]/norm,0*y[1:]]))/E0))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument('eos');p.add_argument('spectra')
    p.add_argument('--step',type=int,default=4);p.add_argument('--height',type=float,default=1.5e-9)
    p.add_argument('--dt',type=float,default=2e-10);p.add_argument('--method',choices=['euler','heun'],default='heun')
    p.add_argument("--verify",action="store_true")
    p.add_argument("--angle",type=int,default=4)
    p.add_argument("--order",type=int,default=5)
    p.add_argument("--cap",type=float,default=0.)
    a=p.parse_args();eos=StandardModel(a.eos)
    row=next(r for r in json.loads(Path(a.spectra).read_text())['rows'] if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
    print(json.dumps(run(row,eos,a.step,a.height,a.dt,a.method,order=a.order,verify=a.verify,angle=a.angle,cap=a.cap),indent=2,allow_nan=False))
