"""Fixed-volume two-state Bose kinetic box with all six pair channels.

q=p/sqrt(2mT_initial), energies q^2+i*d. Equal inertial masses, retain
internal splitting d=gap/T_initial. Time tau=n_initial sigma_LL vT t.
No condensate, gravity, escape or number-changing relativistic processes.
Positive nodal interpolation conserves represented total number/energy.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF
from scipy.optimize import root
from two_state_halo_equilibrium import box_equilibrium,source,pl
from soft_mixing_activation import scales
from late_activation_overlap import critical_density
from resonant_mediator_inverse import M,C,SEC_GYR
try:
    from elastic_acceleration import compiled_rates
except ImportError:
    compiled_rates=None

PAIRS=[(0,0),(0,1),(1,1)]
def geometry(n,qmin=1e-4,qmax=8.,mesh='focused'):
    if mesh=='legacy':
        nl=n//3
        return np.r_[0.,np.geomspace(qmin,.3,nl),np.linspace(.3,qmax,n-nl+1)[1:]]
    nl=max(5,n//6);nt=max(5,n//8);nm=n-nl-nt
    return np.r_[0.,np.geomspace(qmin,.05,nl),np.linspace(.05,4.,nm+1)[1:],np.geomspace(4.,qmax,nt+1)[1:]]

class Box:
    def __init__(self,n,d,u,qmin=1e-4,qmax=8.,alpha=1e-6,mesh='focused'):
        self.e=geometry(n,qmin,qmax,mesh);self.mesh=mesh;l,h=self.e[:-1],self.e[1:]
        self.v=(h**3-l**3)/(6*np.pi**2)
        self.k=.6*(h**5-l**5)/(h**3-l**3);q=np.sqrt(self.k)
        self.n=len(q);n=self.n
        self.d=d;self.u=u;self.alpha=alpha
        self.vol=np.tile(self.v,2);self.energy=np.r_[self.k,self.k+d]
        self.f0=np.r_[1/np.expm1(self.k+alpha),np.zeros(n)]
        self.N0=float(self.vol@self.f0);self.E0=float((self.vol*self.energy)@self.f0)
        dk=np.diff(self.e**2);rows=[];coefs=[];rates=[];labels=[]
        # Integrate all outgoing q3 energies. The incoming triangle for
        # identical particles is doubled off-diagonal. Cross section for
        # distinct final states is twice identical final at fixed amplitude.
        # Half factor for same channel eliminates forward/reverse duplication.
        for a,(s1,s2) in enumerate(PAIRS):
            for b,(s3,s4) in enumerate(PAIRS[a:],a):
                delta=(s1+s2-s3-s4)*d
                amp=((1-u)/u)**(s1+s2+s3+s4)
                pref=amp*(2/(1+(s3==s4)))/(1+(s1==s2))/(2 if a==b else 1)
                for i in range(n):
                    for j in range(i if s1==s2 else 0,n):
                        S=self.k[i]+self.k[j]+delta
                        lo=max(self.k[0],S-self.k[-1]);hi=min(self.k[-1],S-self.k[0])
                        if hi<=lo:continue
                        L=np.maximum(self.k[:-1],lo);H=np.minimum(self.k[1:],hi)
                        wid=np.maximum(H-L,0.);c=(L+H)/2
                        frac=np.clip((c-self.k[:-1])/np.diff(self.k),0.,1.)
                        hats=np.r_[wid*(1-frac),0.]+np.r_[0.,wid*frac]
                        k=np.flatnonzero(hats>0);k4=S-self.k[k]
                        valid=(k4>=self.k[0])&(k4<=self.k[-1]);k=k[valid];k4=k4[valid]
                        h4=np.clip(np.searchsorted(self.k,k4)-1,0,n-2)
                        w=(k4-self.k[h4])/(self.k[h4+1]-self.k[h4])
                        angle=np.maximum(np.minimum(q[i]+q[j],q[k]+np.sqrt(k4))-np.maximum(abs(q[i]-q[j]),abs(q[k]-np.sqrt(k4))),0.)
                        rate=pref*dk[i]*dk[j]*hats[k]*angle*(2 if s1==s2 and i!=j else 1)/(32*np.pi**4*self.N0)
                        ok=rate>0;k,h4,w,rate=k[ok],h4[ok],w[ok],rate[ok]
                        m=len(k)
                        if m:
                            rows.append(np.array([np.full(m,i+s1*n),np.full(m,j+s2*n),k+s3*n,h4+s4*n,h4+1+s4*n],dtype=np.int32))
                            coefs.append(np.array([np.ones(m)]*3+[1-w,w]));rates.append(rate);labels.append(np.full(m,a*3+b))
        self.idx=np.concatenate(rows,axis=1);self.coef=np.concatenate(coefs,axis=1);self.rate=np.concatenate(rates)
        self.labels=np.concatenate(labels);self.sto=self.coef.copy();self.sto[:2]*=-1
        self.inv=1/self.vol
        self.event_energy_error=float(max(abs(np.sum(self.sto*self.energy[self.idx],axis=0))))
    def f(self,y,jac=False,compiled=True):
        if compiled and compiled_rates is not None:
            r,J=compiled_rates(y,self.idx,self.coef,self.sto,self.rate,self.inv,1e-60,jac)
            return J if jac else r
        ff=np.maximum(y,1e-60);g=ff[self.idx];lf=np.log(g);r=np.log1p(1/g)
        rin=np.sum(self.coef[:2]*r[:2],axis=0);rout=np.sum(self.coef[2:]*r[2:],axis=0)
        prod=np.sum(self.coef*lf,axis=0)
        diff=rout-rin
        F=self.rate*np.sign(diff)*np.exp(prod+np.maximum(rin,rout))*(-np.expm1(-abs(diff)))
        if not jac:
            return np.bincount(self.idx.ravel(),weights=(self.sto*F).ravel(),minlength=len(y))*self.inv
        fw=self.rate*np.exp(prod+rout);bw=self.rate*np.exp(prod+rin)
        D=np.empty_like(g)
        D[:2]=self.coef[:2]*(F/g[:2]+bw/(g[:2]*(1+g[:2])))
        D[2:]=self.coef[2:]*(F/g[2:]-fw/(g[2:]*(1+g[2:])))
        D[:2]=np.where(g[:2]<1e-3,self.coef[:2]*(fw/g[:2]-bw/(1+g[:2])),D[:2])
        D[2:]=np.where(g[2:]<1e-3,self.coef[2:]*(fw/(1+g[2:])-bw/g[2:]),D[2:])
        D=np.where(y[self.idx]>1e-60,D,0.)
        J=np.zeros((len(y),len(y)))
        for a in range(5):
            for b in range(5):
                J+=np.bincount(self.idx[a]*len(y)+self.idx[b],weights=self.sto[a]*D[b],minlength=len(y)**2).reshape(J.shape)
        return J*self.inv[:,None]
    def validate(self):
        eq=1/np.expm1(.9*self.energy+.2)
        off=eq*(1+.12*np.sin(np.arange(2*self.n)))
        r=self.f(off);er=self.f(eq)
        v=np.cos(np.arange(2*self.n))*off;eps=1e-6
        fd=(self.f(off+eps*v)-self.f(off-eps*v))/(2*eps)
        heavy=float(self.vol[self.n:]@self.f(self.f0)[self.n:]/self.N0)
        return dict(balance_relative=float(self.vol@abs(er)/(self.vol@abs(r))),
            number_residual=float(abs(self.vol@r)/(self.vol@abs(r))),
            energy_residual=float(abs((self.vol*self.energy)@r)/((self.vol*self.energy)@abs(r))),
            entropy_production=float((self.vol*np.log1p(1/off))@r),
            jacobian_relative=float(np.linalg.norm(self.f(off,True)@v-fd)/np.linalg.norm(fd)),
            event_energy_error=self.event_energy_error,event_count=len(self.rate),
            initial_heavy_gain_per_tau=heavy,
            independent_critical_source=source(self.d,self.u,96,96,alpha=self.alpha)['heavy_initial_number_gain_over_n_sigma_vT'])
    def discrete_eq(self):
        def fun(v):
            t,a=np.exp(v)
            f=1/np.expm1(self.energy/t+a)
            return [self.vol@f/self.N0-1,(self.vol*self.energy)@f/self.E0-1]
        v=root(fun,np.log([.78,.03]))
        if not v.success or max(abs(np.asarray(fun(v.x))))>1e-7:return None,None
        t,a=np.exp(v.x);f=1/np.expm1(self.energy/t+a)
        return f,dict(T_over_initial=float(t),alpha=float(a),heavy_fraction=float(self.v@f[self.n:]/self.N0))
    def run(self,end=1.4,rtol=1e-5,wall=150.):
        started=time.monotonic();check=self.validate();print('check',json.dumps(check),flush=True)
        assert check['balance_relative']<1e-11 and check['number_residual']<1e-12 and check['energy_residual']<1e-12 and check['jacobian_relative']<1e-5 and check['entropy_production']>=0
        f0=np.maximum(self.f0,1e-36);norm=np.maximum(f0,1.)
        sol=BDF(lambda t,y:self.f(y*norm)/norm,0.,f0/norm,end,
            jac=lambda t,y:self.f(y*norm,True)*norm[None,:]/norm[:,None],rtol=rtol,atol=1e-10,max_step=end/50)
        hist=[];steps=0;msg='finished';last=started
        def stats(t,y):
            f=y*norm;s=np.maximum(f,1e-100)
            return dict(tau=float(t),heavy_fraction=float(self.v@f[self.n:]/self.N0),number_error=float(self.vol@f/self.N0-1),
               energy_error=float((self.vol*self.energy)@f/self.E0-1),
               entropy_per_particle=float(self.vol@(np.log1p(s)+s*np.log1p(1/s))/self.N0),
               negative_fraction=float(self.vol@np.maximum(-f,0)/self.N0),
               lowest_L_fraction=float(self.v[0]*f[0]/self.N0),
               f_L_max=float(f[:self.n].max()))
        hist.append(stats(0.,sol.y))
        while sol.status=='running':
            prev=sol.t;yp=sol.y.copy();note=sol.step();steps+=1
            st=stats(sol.t,sol.y)
            if st['negative_fraction']>1e-8:
                sol.t=prev;sol.y=yp;msg='negative accepted number';break
            if sol.status=='failed':msg=str(note);break
            if steps%5==0 or sol.status=='finished':hist.append(st)
            if time.monotonic()-last>12:
                print('progress',json.dumps(st),flush=True);last=time.monotonic()
            if time.monotonic()-started>wall:msg='wall budget';break
        if hist[-1]['tau']!=sol.t:hist.append(stats(sol.t,sol.y))
        f_eq,eq=self.discrete_eq()
        if eq is not None: eq['final_number_L1_relative']=float(self.vol@abs(sol.y*norm-f_eq)/self.N0)
        return dict(discrete_equilibrium=eq,mesh=self.mesh,success=sol.status=='finished' and msg=='finished',message=msg,n=self.n,d=self.d,u=self.u,alpha=self.alpha,
            N0=self.N0,E0=self.E0,expected_critical_box=box_equilibrium(self.d),
            requested_tau=end,achieved_tau=float(sol.t),checks=check,steps=steps,wall_seconds=time.monotonic()-started,
            final=hist[-1],max_N_error=max(abs(x['number_error']) for x in hist),max_E_error=max(abs(x['energy_error']) for x in hist),
            history=hist,edges=self.e.tolist(),initial=self.f0.tolist(),final_f=(sol.y*norm).tolist(),
            scope='Homogeneous two-state number-conserving ideal NR kinetic box. No condensation, halo transport, finite-volume mass shift or 4->2.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--n',type=int,default=48);p.add_argument('--V0',type=float,default=44.1)
    p.add_argument('--mesh',choices=['focused','legacy'],default='focused');p.add_argument('--sigma',type=float,default=.3);p.add_argument('--tau',type=float);p.add_argument('--alpha',type=float,default=1e-6);p.add_argument('--rtol',type=float,default=1e-5);p.add_argument('--wall',type=float,default=120.);a=p.parse_args()
    par=scales(a.sigma,z_on=8.);vT=np.sqrt(2)*a.V0*1e5;T=M*(vT/C)**2/2
    tau=a.tau or critical_density(vT)*a.sigma*vT*10*SEC_GYR
    box=Box(a.n,par['mass_gap_eV']/T,par['u'],alpha=a.alpha,mesh=a.mesh)
    r=box.run(tau,a.rtol,a.wall);r.update(params=par,V0=a.V0,T_eV=T)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k not in ['history','initial','edges','final_f','params']}),flush=True)
