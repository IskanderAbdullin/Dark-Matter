"""Closed energy/count moment comparator for a shared thermal photon bath."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root
from scipy.special import expit, roots_legendre

from early_pole_hard_gauge_pair_source_reference import calc
from early_pole_joint_fermion_thermal_reference import Joint


class Pair:
    def __init__(self,mass,radial,incident,final,azimuth):
        q,w=roots_legendre(radial)
        p=.5+19.75*(q+1)
        G=19.75*w*p*p/(np.pi*np.pi)
        p1,p2=np.meshgrid(p,p,indexing='ij')
        vol=G[:,None]*G[None,:]
        z,v=roots_legendre(incident)
        ss0=4*p1*p2
        a=np.minimum(1.,np.maximum(1.,4*mass*mass)/ss0)
        u=a[:,:,None]+(1-a[:,:,None])*(1+z)/2
        w0=vol[:,:,None]*(1-a[:,:,None])*v*u
        ss=ss0[:,:,None]*u
        total=(p1+p2)[:,:,None]+np.zeros_like(ss)
        boost=np.sqrt(np.maximum(total*total-ss,0.))
        alpha=np.divide((p1-p2)[:,:,None],boost,out=np.zeros_like(ss),where=boost>0.)
        alpha=np.clip(alpha,-1.,1.)
        r=np.minimum(1.,4*mass*mass/ss)
        beta=np.sqrt(1-r)
        th,tw=roots_legendre(final)
        az,aw=roots_legendre(azimuth)
        th*=.8;tw*=.8
        accum=[]
        for cos,wcos in zip(th,tw,strict=True):
            v0=beta*cos
            den=1-v0*v0
            dcs=beta*(4*(1+r)/den-2-4*r*r/(den*den))/(16*np.pi*ss)
            for phi,wphi in zip(np.pi*(az+1),aw/2,strict=True):
                ca=alpha*cos+np.sqrt((1-alpha**2)*(1-cos*cos))*np.cos(phi)
                E0=total/2+boost*beta*ca/2
                E1=total-E0
                K=w0*dcs*wcos*wphi
                keep=K>1e-35
                assert min(E0[keep].min(initial=mass),E1[keep].min(initial=mass))>=mass-1e-10
                accum.append(np.array([np.broadcast_to(p1[:,:,None],ss.shape)[keep],
                                       np.broadcast_to(p2[:,:,None],ss.shape)[keep],E0[keep],E1[keep],K[keep]]))
        self.g1,self.g2,self.f1,self.f2,self.w=np.concatenate(accum,axis=1)
        self.E=self.g1+self.g2
        assert max(abs(self.f1+self.f2-self.E))<1e-13
        self.bare=self.w@((1/np.expm1(self.g1))*(1/np.expm1(self.g2)))

    def rate(self,gamma_temperature,fermion_temperature,mu_over_fermion_T):
        # Factor the reverse term using the multiplier from detailed balance.
        th=gamma_temperature
        tf=fermion_temperature
        a=mu_over_fermion_T
        b1=self.g1/th
        b2=self.g2/th
        with np.errstate(divide='raise'):
            log_bose1=-(b1+np.log(-np.expm1(-b1)))
            log_bose2=-(b2+np.log(-np.expm1(-b2)))
        forward=(np.log(self.w)+log_bose1+log_bose2-
                 np.logaddexp(0.,a-self.f1/tf)-np.logaddexp(0.,a-self.f2/tf))
        L=2*a+self.E*(1/th-1/tf)
        net=np.sign(-L)*np.exp(forward+np.maximum(L,0.))*(-np.expm1(-abs(L)))
        return np.array([net.sum(),.5*net@self.E])

    def heat_control(self):
        out=[]
        for Tg,Tf in ((1.,.8),(.85,1.05)):
            chemical=brentq(lambda a:self.rate(Tg,Tf,a)[0],-6.,6.,xtol=3e-14)
            n,heat=self.rate(Tg,Tf,chemical)
            assert abs(n)<1e-12 and heat*(Tg-Tf)>0
            entropy=(1./Tf-1./Tg)*heat-chemical*n
            assert entropy>0
            out.append(dict(photon_T=Tg,fermion_T=Tf,mu_over_fermion_T=chemical,
                            count_residual=n,heat_at_zero_count=heat,
                            two_bath_entropy_production=entropy))
        return out


class Model:
    def __init__(self,s,theta,g,order,fixed):
        x=s['x_end'];rb=s['mass_ratio']
        self.fixed=fixed;self.g=g;self.x=x
        self.theta0=theta
        self.H0=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
        self.photp=np.pi**2/15
        self.Phot0=self.photp*theta**4
        self.masses=(rb*x,x)
        self.pairs=[Pair(m,*order) for m in self.masses]
        for A,m in zip(self.pairs,self.masses,strict=True):
            assert abs(A.bare/calc(m,1.,.8,80)-1)<.008
        G,W=roots_legendre(96)
        self.q=20*(1+G)
        e=np.concatenate([np.hypot(self.q,m) for m in self.masses])
        w=np.tile(4*20*W*self.q**2/(2*np.pi**2),2)
        u=np.arange(len(e))<len(self.q)
        self.features=np.column_stack((u,~u,e))
        self.B=(w[:,None]*self.features).T
        self.ref_density=np.pi**2*106.75/30
        shapes={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
        df=Joint(shapes[tuple(s['cells'])],rb,0.,'none')
        fd=np.array(s['endpoint_fermion_occupations'])
        fp=np.array(s['endpoint_parent_occupations'])
        self.y0=np.array([df.d.wf@fd,df.pw@fp,df.d.wf@(df.daughter_energies(x)*fd)+
                                               df.pw@(df.parent_energies(x)*fp)])
        self.E0=self.y0[2]+self.Phot0

    def targets(self,y):
        scale=np.maximum(abs(y),1e-12)
        F=self.features;B=self.B
        z0=[np.log(max(y[i],1e-20)/.365) for i in range(2)]+[-1.]
        def fd(z):
            return expit(F@z)
        def derivative(z):
            q=F@z
            return F*(expit(q)*expit(-q))[:,None]
        sol=root(lambda z:(B@fd(z)-y)/scale,z0,
                 jac=lambda z:(B@derivative(z))/scale[:,None],tol=2e-11)
        assert np.max(abs((B@fd(sol.x)-y)/scale))<2e-11,(sol.message,y,sol.x)
        return sol.x[:2],-1/sol.x[-1]

    def photon_temperature(self,y):
        return self.theta0 if self.fixed else ((self.E0-y[2])/self.photp)**.25

    def rhs(self,_,y):
        f,tf=self.targets(y)
        th=self.photon_temperature(y)
        a,b=[q.rate(th,tf,mu) for q,mu in zip(self.pairs,f,strict=True)]
        return self.g**4/self.H0*np.array([a[0],b[0],a[1]+b[1]])

    def equilibrium(self):
        def moments(th):
            return self.B@expit(-self.features[:,2]/th)
        if self.fixed:
            T=self.theta0
        else:
            T=brentq(lambda q:moments(q)[2]+self.photp*q**4-self.E0,.01,4.,xtol=1e-14)
        a=moments(T)
        return dict(T=T,daughter_number=a[0],parent_number=a[1],fermion_energy=a[2],
                    photon_energy=self.photp*T**4)

    def go(self,end):
        ts=np.linspace(0,end,25)
        sol=solve_ivp(self.rhs,(0,end),self.y0,method='DOP853',t_eval=ts,rtol=3e-8,atol=3e-11)
        assert sol.success
        result=[]
        for t,y in zip(sol.t,sol.y.T):
            a,tf=self.targets(y)
            th=self.photon_temperature(y)
            photon=self.photp*th**4
            rates=self.rhs(t,y)
            entropy=(1/tf-1/th)*rates[2]-a@rates[:2]
            assert entropy>=-3e-13
            result.append(dict(reference_time=t,number_daughter=y[0],number_parent=y[1],
                               fermion_energy=y[2],photon_energy=photon,photon_T=th,
                               fermion_T=tf,chemical_difference=(a[1]-a[0])*tf,
                               energy_defect=photon+y[2]-self.E0,paired_entropy_production=entropy))
        if not self.fixed:assert max(abs(v['energy_defect']) for v in result)<3e-14
        for A in self.pairs:
            assert np.max(abs(A.rate(.93,.93,0.)))==0.
        return dict(x=self.x,g=self.g,fixed_temperature_photons=self.fixed,
                    unblocked_vacuum_hard_moment_check=[A.bare for A in self.pairs],
                    initial_rhs=self.rhs(0,self.y0).tolist(),
                    zero_count_nonzero_heat_control=self.pairs[0].heat_control(),
                    analytic_stationary_moments=self.equilibrium(),
                    rows=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('spec')
    p.add_argument('--theta',type=float,default=1.)
    p.add_argument('--g',type=float,default=.001)
    p.add_argument('--end',type=float,default=.2)
    p.add_argument('--order',type=int,default=14)
    p.add_argument('--refined-angle',action='store_true')
    p.add_argument('--fixed',action='store_true')
    a=p.parse_args()
    f=Path(a.spec)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    s=json.loads(raw)
    d=Model(s,a.theta,a.g,(a.order,14,8,8) if a.refined_angle else (a.order,10,6,6),a.fixed)
    print(json.dumps(d.go(a.end),indent=2,allow_nan=False))
