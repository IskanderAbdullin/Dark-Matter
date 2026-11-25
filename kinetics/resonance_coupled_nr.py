"""Conservative expanding NR A+A <-> R resonance kinetic pilot.

All variables use momentum in units of the original comoving peak.
R rest excess corresponds to 2*b(a)^2 in A kinetic-energy units.
Fractional stoichiometric interpolation conserves A+2R and instantaneous
kinetic+excess energy. Positive geometric factors enforce BE balance.
Truncated cold-sector diagnostic; no relativistic tails or dressing.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF,Radau
from scipy.special import roots_legendre
from quartic_consistency import shell_nodes
from resonant_mediator_inverse import width,averages,HBAR,C,M
from late_activation_overlap import H

def grid(e):
    l,u=e[:-1],e[1:]
    vol=(u**3-l**3)/(6*np.pi**2)
    y=.6*(u**5-l**5)/(u**3-l**3)
    return vol,y

def map_A(row,n,qmax):
    p0=row['final']['peak_bin_q']
    e=np.asarray(row['frozen_cell_spectrum']['edges'])/p0
    f=np.maximum(row['frozen_cell_spectrum']['f_phi'],0.)
    # Include zero and exact upper face; preserve input shell number when merging.
    sub=np.unique(np.r_[0.,e[e<qmax],qmax])
    l,u=sub[:-1],sub[1:];mid=(l+u)/2
    j=np.searchsorted(e,mid,side='right')-1
    ff=np.where((j>=0)&(j<len(f)),f[np.clip(j,0,len(f)-1)],0.)
    v,y=grid(sub);N=ff*v
    ix=np.unique(np.rint(np.linspace(0,len(ff),min(n,len(ff))+1)).astype(int))
    e2=sub[ix];v2,y2=grid(e2)
    g=np.array([N[a:b].sum() for a,b in zip(ix[:-1],ix[1:])])/v2
    fullN=float(shell_nodes(e*p0,f,12)[1].sum()/p0**3)
    return e2,g,dict(full_number=fullN,retained_fraction=float(v2@g/fullN),
                       remesh_energy_relative_error=float((v2*y2)@g/((v*y)@ff)-1))

class Model:
    def __init__(self,eA,eR,vr,v0,gamma,order=2,regularization=1e-12):
        self.eA=eA;self.eR=eR
        self.va,self.ya=grid(eA);self.vr,self.yr=grid(eR)
        self.na=len(self.va);self.nr=len(self.vr);self.n=self.na+self.nr
        self.qr=np.sqrt(self.yr);self.yr*=.5
        self.vol=np.r_[self.va,self.vr];self.inv=1/self.vol
        self.charge=self.vol*np.r_[np.ones(self.na),2*np.ones(self.nr)]
        self.vres=vr;self.v0=v0;self.gamma=gamma;self.order=order;self.key=None
        self.u,self.w=roots_legendre(order);self.floor=1e-90;self.regularization=regularization
    def geometry(self,z):
        if self.key==z:return
        self.key=z;self.b=self.vres/(2*self.v0*np.exp(-z));b=self.b
        a=self.qr**2/4+b*b;d=self.qr*b
        lim=np.clip(np.minimum((self.ya[-1]-a)/d,(a-self.ya[0])/d),0.,1.)
        bounds=np.sort(np.clip(np.c_[np.zeros(self.nr),lim,
            (self.ya[None,:]-a[:,None])/d[:,None],(a[:,None]-self.ya[None,:])/d[:,None]],0.,lim[:,None]),axis=1)
        lo,hi=bounds[:,:-1],bounds[:,1:];width=hi-lo
        jj,kk=np.nonzero(width>0)
        if not len(jj):
            self.idx=np.empty((5,0),dtype=np.int32);self.c=np.empty((5,0));self.s=self.c.copy();self.rate=np.empty(0);self.accepted=lim
            self.energy_error=0.;return
        W=width[jj,kk];mid=(lo[jj,kk]+hi[jj,kk])/2
        c=mid[:,None]+W[:,None]*self.u/2
        R=np.repeat(jj,self.order);c=c.ravel()
        y1=a[R]+d[R]*c;y2=a[R]-d[R]*c
        def interp(y):
            i=np.clip(np.searchsorted(self.ya,y)-1,0,self.na-2)
            w=np.clip((y-self.ya[i])/(self.ya[i+1]-self.ya[i]),0.,1.)
            return i,w
        i,wi=interp(y1);j,wj=interp(y2)
        self.idx=np.array([self.na+R,i,i+1,j,j+1],dtype=np.int32)
        self.c=np.array([np.ones(len(R)),1-wi,wi,1-wj,wj])
        self.s=self.c.copy();self.s[0]*=-1
        self.rate=(self.vr[R]*np.repeat(W,self.order)*np.tile(self.w,len(jj))/2)*self.gamma/H(np.exp(-z)-1)
        self.accepted=lim
        self.energy_error=float(np.max(abs((1-wi)*self.ya[i]+wi*self.ya[i+1]+(1-wj)*self.ya[j]+wj*self.ya[j+1]-(self.yr[R]+2*b*b))))
    def flux(self,z,f,jac=False):
        self.geometry(z);ff=np.maximum(f,self.floor);lf=np.log(ff);r=np.log1p(1/ff)
        vals=ff[self.idx];logprod=np.sum(self.c*lf[self.idx],axis=0)
        rin=r[self.idx[0]];rout=np.sum(self.c[1:]*r[self.idx[1:]],axis=0)
        delta=rout-rin
        F=self.rate*np.sign(delta)*np.exp(logprod+np.maximum(rin,rout))*(-np.expm1(-abs(delta)))
        # Fractional powers at an empty node otherwise have a singular
        # derivative. A C1 positive continuation below a tiny occupation
        # gives exact inward boundary flux. Use stable Bose subtraction
        # on all ordinary events, and explicit products only on the
        # near-empty subset.
        low=np.any((f[self.idx]<self.regularization)&(self.c>0),axis=0) if self.regularization else np.zeros(len(F),bool)
        if jac:
            fw=self.rate*np.exp(logprod+rout);bw=self.rate*np.exp(logprod+rin)
            D=np.empty_like(self.c)
            D[0]=F/vals[0]+bw/(vals[0]*(1+vals[0]))
            D[1:]=self.c[1:]*(F[None,:]/vals[1:]-fw[None,:]/(vals[1:]*(1+vals[1:])))
            near=vals<1e-3
            D[0]=np.where(near[0],fw/vals[0]-bw/(1+vals[0]),D[0])
            D[1:]=np.where(near[1:],self.c[1:]*(fw[None,:]/(1+vals[1:])-bw[None,:]/vals[1:]),D[1:])
            D=np.where(f[self.idx]>self.floor,D,0.)
        if np.any(low):
            c=self.c[:,low];raw=f[self.idx[:,low]];g=np.maximum(raw,0.);fl=self.regularization
            gsafe=np.maximum(g,fl)
            av=np.exp(c*np.log(gsafe));da=av*c/gsafe
            t=np.clip(g/fl,0.,1.);mask=(g<fl)&(c>0)
            scale=fl**c
            av=np.where(mask,scale*((2-c)*t-(1-c)*t*t),av)
            da=np.where(mask,scale*((2-c)-2*(1-c)*t)/fl,da)
            bl=np.exp(c*np.log1p(g));db=bl*c/(1+g)
            da=np.where(raw>=0,da,0.);db=np.where(raw>=0,db,0.)
            A=np.vstack((av[0],bl[1:]));B=np.vstack((bl[0],av[1:]))
            rate=self.rate[low]
            F[low]=rate*(np.prod(A,axis=0)-np.prod(B,axis=0))
            if jac:
                Ap=np.vstack((da[0],db[1:]));Bp=np.vstack((db[0],da[1:]))
                for k in range(5):
                    ix=np.arange(5)!=k
                    D[k,low]=rate*(Ap[k]*np.prod(A[ix],axis=0)-Bp[k]*np.prod(B[ix],axis=0))
        if not jac:return np.bincount(self.idx.ravel(),weights=(self.s*F).ravel(),minlength=self.n)*self.inv
        J=np.zeros((self.n,self.n))
        for a in range(5):
            for b in range(5):
                J+=np.bincount(self.idx[a]*self.n+self.idx[b],weights=self.s[a]*D[b],minlength=self.n*self.n).reshape(self.n,self.n)
        return J*self.inv[:,None]
    def energy(self,z):
        self.geometry(z);return self.vol*np.r_[self.ya,self.yr+2*self.b*self.b]
    def validate(self,z,fA):
        self.geometry(z)
        eq=1/np.expm1(.002*np.r_[self.ya,self.yr+2*self.b*self.b]+np.r_[np.full(self.na,.7),np.full(self.nr,1.4)])
        d=self.flux(z,eq)
        test=eq*(1+.25*np.sin(np.arange(self.n)))
        d2=self.flux(z,test);J=self.flux(z,test,True)
        v=np.cos(np.arange(self.n))*test;eps=1e-6
        fd=(self.flux(z,test+eps*v)-self.flux(z,test-eps*v))/(2*eps)
        a=np.r_[fA,np.zeros(self.nr)]
        inv=self.flux(z,a)[self.na:]
        F,G=averages(self.qr,self.b,self.eA**2,fA)
        ref=self.gamma/H(np.exp(-z)-1)*F
        norm=max(self.vr@ref,1e-300)
        return dict(balance_relative=float(self.vol@abs(d)/max(self.vol@abs(d2),1e-300)),
            number_residual=float(abs(self.charge@d2)/max(self.charge@abs(d2),1e-300)),
            energy_residual=float(abs(self.energy(z)@d2)/max(self.energy(z)@abs(d2),1e-300)),
            jacobian_relative=float(np.linalg.norm(J@v-fd)/max(np.linalg.norm(fd),1e-300)),
            entropy_production=float((self.vol*np.log1p(1/test))@d2),
            direct_inverse_relative_error=float((self.vr@inv)/norm-1),event_energy_error=self.energy_error)

def run(row,vr_kms,nA=100,nR=80,qmax=20.,rstart=2e6,bend=10.,rtol=3e-5,order=2,wall=160.,regularization=1e-12,method='BDF'):
    start=time.monotonic();h=row['final'];p0=h['peak_bin_q']
    v0=p0*2.7255*8.617333262e-5*(3.909/row['g_bath'])**(1/3)/h['theta']/M*C
    params=width(vr_kms);gamma=params['elastic_width_eV']/HBAR
    rmax=min(rstart,.03*C/(qmax*v0))
    rend=vr_kms*1e5/(2*v0*bend)
    z0,z1=-np.log(rmax),-np.log(rend)
    assert z1>z0
    eA,fA,meta=map_A(row,nA,qmax)
    eR=np.r_[0.,np.geomspace(1e-4,2*qmax,nR+1)]
    m=Model(eA,eR,vr_kms*1e5,v0,gamma,order,regularization)
    N0=float(m.va@fA);E0=float((m.va*m.ya)@fA)
    check=m.validate(z0,fA)
    assert check['balance_relative']<1e-10 and check['jacobian_relative']<1e-5 and check['number_residual']<1e-12 and check['energy_residual']<1e-12 and check['entropy_production']>=0,check
    f0=np.r_[fA,np.zeros(m.nr)]
    # A tiny weighted scale prevents empty soft cells from dominating absolute error.
    norm=np.maximum(np.maximum(f0,1.),N0/(m.n*m.vol)*1e-5)
    y0=np.r_[np.maximum(f0,1e-38)/norm,0.]
    def rhs(z,y):
        f=y[:-1]*norm;d=m.flux(z,f)
        return np.r_[d/norm,4*m.b*m.b*(m.vr@f[m.na:])/E0]
    def jac(z,y):
        f=y[:-1]*norm;J=m.flux(z,f,True)*norm[None,:]/norm[:,None]
        out=np.zeros((m.n+1,m.n+1));out[:-1,:-1]=J
        out[-1,m.na:-1]=4*m.b*m.b*m.vr*norm[m.na:]/E0
        return out
    cls={"BDF":BDF,"Radau":Radau}[method]
    # Shift integration time to start at zero; geometry still uses ln(a).
    zzrhs=lambda t,y:rhs(t+z0,y)
    zzjac=lambda t,y:jac(t+z0,y)
    sol=cls(zzrhs,0.,y0,z1-z0,jac=zzjac,rtol=rtol,atol=np.r_[np.full(m.n,1e-10),1e-12],max_step=.025)
    hist=[];steps=0;msg='finished';last=start
    def stats(z,y):
        f=y[:-1]*norm;A,R=f[:m.na],f[m.na:];m.geometry(z)
        s=np.maximum(f,1e-300)
        return dict(z=float(np.exp(-z)-1),b=float(m.b),R_charge_fraction=float(2*m.vr@R/N0),
            charge_relative_error=float(m.charge@f/N0-1),
            energy_continuity_residual=float(m.energy(z)@f/E0-1-y[-1]),
            entropy_per_retained_charge=float(m.vol@(np.log1p(s)+s*np.log1p(1/s))/N0),
            A_q2_per_initial=float((m.va*m.ya)@A/E0),
            A_hot_fraction_over_initial=float(m.va[m.ya>36]@A[m.ya>36]/N0),
            max_fA=float(A.max()),max_fR=float(R.max()),
            R_dropped_angle_weighted_fraction=float(m.vr@(R*(1-m.accepted))/max(m.vr@R,1e-300)),
            negative_charge_fraction=float(m.charge@np.maximum(-f,0)/N0))
    hist.append(stats(sol.t+z0,sol.y))
    while sol.status=='running':
        zp=sol.t;yp=sol.y.copy();note=sol.step();steps+=1
        f=sol.y[:-1]*norm
        if m.charge@np.maximum(-f,0)/N0>1e-8:
            sol.t=zp;sol.y=yp;msg='negative accepted charge';break
        if sol.status=='failed':msg=str(note);break
        if steps%10==0 or sol.status=='finished':
            hist.append(stats(sol.t+z0,sol.y))
            if time.monotonic()-last>8 or sol.status=='finished':
                print('progress',json.dumps(hist[-1]),flush=True);last=time.monotonic()
        if time.monotonic()-start>wall:msg='wall budget';break
    if hist[-1]['z']!=float(np.exp(-(sol.t+z0))-1):hist.append(stats(sol.t+z0,sol.y))
    f=sol.y[:-1]*norm
    return dict(success=sol.status=='finished' and msg=='finished',message=msg,vr_kms=vr_kms,
        nA=m.na,nR=m.nr,qmax=qmax,rstart=rmax,rend=rend,rtol=rtol,order=order,params=params,
        regularization=regularization,method=method,meta=meta,checks=check,steps=steps,nfev=sol.nfev,njev=sol.njev,wall_seconds=time.monotonic()-start,
        final=hist[-1],maximum_R_fraction=max(v['R_charge_fraction'] for v in hist),
        maximum_charge_residual=max(abs(v['charge_relative_error']) for v in hist),
        maximum_energy_residual=max(abs(v['energy_continuity_residual']) for v in hist),
        maximum_dropped_R_decay_fraction=max(v['R_dropped_angle_weighted_fraction'] for v in hist),
        history=hist,eA=eA.tolist(),eR=eR.tolist(),initial_fA=fA.tolist(),final_fA=f[:m.na].tolist(),final_fR=f[m.na:].tolist(),
        scope='Cold truncated NR A+ A <-> R, no hot interactions, vacuum width, represented FRW energy. Not a full relic prediction.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--vr',type=float,default=90.);p.add_argument('--nA',type=int,default=100);p.add_argument('--nR',type=int,default=80)
    p.add_argument('--qmax',type=float,default=20.);p.add_argument('--rstart',type=float,default=2e6)
    p.add_argument('--bend',type=float,default=10.);p.add_argument('--rtol',type=float,default=3e-5)
    p.add_argument('--method',choices=['BDF','Radau'],default='BDF');p.add_argument('--regularization',type=float,default=1e-12);p.add_argument('--order',type=int,default=2);p.add_argument('--wall',type=float,default=160.);a=p.parse_args()
    r=run(json.loads(Path(a.input).read_text())[0],a.vr,a.nA,a.nR,a.qmax,a.rstart,a.bend,a.rtol,a.order,a.wall,a.regularization,a.method)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print('FINAL',json.dumps({k:v for k,v in r.items() if k not in ('history','eA','eR','initial_fA','final_fA','final_fR')}),flush=True)
