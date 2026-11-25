"""Coupled portal production with state-dependent Hartree quasiparticle mass.

Uses U=integral E(m_eff) f - lambda V^2/8 and trace=m_vac^2 V. The gap is
stationary, so dU/df is exactly the shell quasiparticle energy. Each decay
and optional elastic event conserves this derivative at the current state.
For Newton iteration freeze the gap in the collision graph derivative;
the actual RHS recomputes it. No zero mode, coherent field, 2<->4, thermal
vertex correction or proof of bare kinetic validity. Stop at lambda*f guard.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.integrate import BDF
from scipy.optimize import brentq
from scipy.special import expit
from parent_angle_kernel import ParentAngle
from laine_portal_rates import PortalTable,make_models
from elastic_implicit_trial import ImplicitElastic
from quartic_consistency import lambda_from_sigma,MPL
from radiation_fixed_bath import omega_from_number
from two_mass_bath import invariants

_,TableResolved=make_models()

class HartreePortal(ParentAngle,TableResolved):
    def __init__(self,sigma=0.,elastic=False,**kw):
        super().__init__(**kw)
        self.sigma=sigma;self.lam=float(lambda_from_sigma(sigma));self.do_elastic=elastic and sigma>0
        self.mass_key=None;self.graph=None
        self.ec=self.lam**2/(2048*np.pi**5)*MPL/(1.66*np.sqrt(self.g_bath)*self.parent_GeV)
    def gap_state(self,x,p):
        bare=x*self.scalar_ratio
        pp=np.maximum(p,0.)*self.vol
        def var(m):return float(pp@np.sum(self.ws/np.hypot(self.ps,m),axis=1))
        V0=var(bare);delta=self.lam*V0/2
        if delta>0:
            r=brentq(lambda u:u-var(np.sqrt(bare*bare+u*delta))/V0,0.,1+1e-12,xtol=1e-13)
            m=np.sqrt(bare*bare+r*delta)
        else:m=bare
        V=var(m);return m,V,self.lam*V*V/8
    def prepare(self,x,p):
        m,V,B=self.gap_state(x,p)
        if m+x*self.fermion_ratio>=x:raise ValueError('Hartree gap closed decay channel; outside pilot')
        if self.mass_key!=(x,m):
            bare=self.scalar_ratio;kappa=self.kappa
            _,D=invariants(x,m,x*self.fermion_ratio)
            _,D0=invariants(x,bare*x,x*self.fermion_ratio)
            # Hold bare spin-summed matrix element fixed; phase space ratio only.
            self.kappa=kappa*D/D0
            self.scalar_ratio=m/x;self.geom_x=-1.
            try:ParentAngle.geometry(self,x)
            finally:self.scalar_ratio=bare;self.kappa=kappa
            self.traces[:self.n]=self.vol*np.sum(self.ws*(bare*x)**2/np.hypot(self.ps,m),axis=1)
            self.trace_weights[:self.n]=self.traces[:self.n]
            if self.do_elastic and (self.graph is None or self.graph.m!=m):
                self.graph=ImplicitElastic(self.edges,m)
                assert np.max(abs(self.energy_phi-self.graph.K-m)/self.energy_phi)<1e-11
            self.mass_key=(x,m)
        return m,V,B
    def eval_hartree(self,z,y,jac=False):
        x=np.exp(z);yd=y[:self.dim];theta=np.exp(y[-4]);rb=self.rad0*theta**4
        m,V,B=self.prepare(x,yd[:self.n])
        total=rb+self.energy_weights@yd-B;Hf=np.sqrt(self.rad0/total)
        rate,rateT=self.rates(x,theta,Hf);eq=expit(-self.energies/theta);bath=rate*(eq-yd)
        if jac:dark,J=self.dark(yd,Hf,True)
        else:dark=self.dark(yd,Hf)
        if self.do_elastic:
            c=self.ec*x*Hf
            dark[:self.n]+=c*self.graph.rates_jac(yd[:self.n],False)
            if jac:J[:self.n,:self.n]+=c*self.graph.rates_jac(yd[:self.n],True)
        dy=dark+bath;dt=-self.energy_weights@bath/(4*rb)
        if not jac:return np.r_[dy,dt,self.traces@yd,self.source1@bath,self.source2@bath]
        # Approximate Newton matrix: exact frozen-mass collision derivatives,
        # plus exact stationary-energy derivative for H and bare trace gradient.
        dH=-.5/total
        J-=np.diag(rate);J+=dH*np.outer(dy,self.energy_weights)
        bathT=rate*eq*(1-eq)*self.energies/theta+rateT*(eq-yd)
        JT=bathT+dH*dy*4*rb
        BB=-np.diag(rate)+dH*np.outer(bath,self.energy_weights)
        BT=bathT+dH*bath*4*rb
        jtemp=-self.energy_weights@BB/(4*rb);tempT=-self.energy_weights@BT/(4*rb)-4*dt
        traces=self.traces.copy()
        V3=(self.vol*np.maximum(yd[:self.n],0.))@np.sum(self.ws/np.hypot(self.ps,m)**3,axis=1)
        traces[:self.n]/=1+self.lam*V3/4
        D=np.column_stack((J,JT,np.zeros((self.dim,3))))
        return np.vstack((D,np.r_[jtemp,tempT,0.,0.,0.],np.r_[traces,0.,0.,0.,0.],
                          np.r_[self.source1@BB,self.source1@BT,0.,0.,0.],np.r_[self.source2@BB,self.source2@BT,0.,0.,0.]))
    def check_local(self,x,y):
        p=y[:self.n];m,V,B=self.prepare(x,p)
        yd=y[:self.dim];theta=np.exp(y[-4]);rb=self.rad0*theta**4
        total=rb+self.energy_weights@yd-B
        # Direct independent derivative of stationary energy functional.
        dr=np.random.default_rng(123).uniform(-1,1,size=self.n)*np.maximum(p,1e-20)
        def U(pp):
            mass,v,b=self.gap_state(x,pp)
            return float((self.vol*np.sum(self.ws*np.hypot(self.ps,mass),axis=1))@pp-b)
        eps=1e-5
        dU=(U(p+eps*dr)-U(p-eps*dr))/(2*eps)
        exact=self.energy_weights[:self.n]@dr
        return dict(energy_gradient_relative_error=float(abs(dU-exact)/max(abs(exact),1e-300)))
    def run(self,xend,rtol,wall,guard):
        start=time.monotonic();x0=.1;y0=np.zeros(self.dim+4)
        atol=np.r_[np.full(self.n,1e-12),np.full(self.n+self.n1,1e-11),1e-13,1e-12,1e-13,1e-13]
        sol=BDF(lambda z,y:self.eval_hartree(z,y),np.log(x0),y0,np.log(xend),
            jac=lambda z,y:self.eval_hartree(z,y,True),rtol=rtol,atol=atol,max_step=.035)
        hist=[];psd=[];steps=0;message='finished';last=start
        targets=iter([x for x in [.125,.15,.175,.2,.3,.5,1.,2.,3.,5.,10.,20.,40.,80.] if x<xend])
        target=next(targets,None)
        def state(z,y):
            x=np.exp(z);yd=y[:self.dim];p=yd[:self.n];g=yd[2*self.n:];f=yd[self.n:2*self.n]
            m,V,B=self.prepare(x,p);theta=np.exp(y[-4]);rb=self.rad0*theta**4
            N=float(self.vol@p);N1=float(self.v1@g);N2=float(self.v2@f)
            U=float(self.energy_weights@yd-B);P=(U-self.traces@yd)/3
            cold=np.maximum(np.minimum(self.edges[1:],.1*theta)**3-self.edges[:-1]**3,0)/(self.edges[1:]**3-self.edges[:-1]**3)
            return dict(x=float(x),number=N,omega_proxy=float(omega_from_number(N/theta**3,gs=self.g_bath)),
                m_eff=m,lambda_fmax=float(self.lam*p.max()),fraction_lambda_f_gt1=float(self.vol[p*self.lam>1]@p[p*self.lam>1]/max(N,1e-300)),
                low_two_fraction=float(self.vol[:2]@p[:2]/max(N,1e-300)),hot_fraction=float((p*self.vol)@(1-cold)/max(N,1e-300)),
                theta=float(theta),mean_field_B=B,phi_energy=float(self.energy_weights[:self.n]@p-B),
                number_parent_residual=float((N+N1-y[-2])/max(N,1e-20)),
                number_daughter_residual=float((N2-N-y[-1])/max(N,1e-20)),
                continuity_residual=float((rb+U-self.rad0-y[-3])/self.rad0),
                min_occupation=float(yd.min()),max_fermion=float(yd[self.n:].max()),elapsed=time.monotonic()-start)
        while sol.status=='running':
            msg=sol.step();steps+=1
            if sol.status=='failed':message=str(msg);break
            x=np.exp(sol.t);y=sol.y
            if np.min(y[:self.dim])<-1e-7 or np.max(y[self.n:self.dim])>1+1e-7:message='invalid accepted occupation';break
            stop=self.lam*np.max(y[:self.n])>=guard
            if target is not None and x>=target:
                interp=sol.dense_output()
                while target is not None and x>=target:
                    z=np.log(target);yy=interp(z);h=state(z,yy);hist.append(h);psd.append(yy[:self.n].tolist());print('state',json.dumps(h),flush=True)
                    target=next(targets,None)
            if stop:message='tree kinetic occupation guard';break
            if time.monotonic()-start>wall:message='wall budget';break
            if time.monotonic()-last>15:
                print('progress',float(x),steps,sol.nfev,float(self.lam*np.max(y[:self.n])),flush=True);last=time.monotonic()
        h=state(sol.t,sol.y);hist.append(h);psd.append(sol.y[:self.n].tolist());check=self.check_local(h['x'],sol.y)
        return dict(success=sol.status=='finished' and message=='finished',message=message,n=self.n,nparent=self.n1,
            sigma=self.sigma,elastic=self.do_elastic,guard=guard,rtol=rtol,nfev=sol.nfev,njev=sol.njev,steps=steps,
            final=h,history=hist,scalar_history=psd,edges=self.edges.tolist(),checks=check,
            max_continuity_residual=max(abs(v['continuity_residual']) for v in hist),
            max_number_residual=max(max(abs(v['number_parent_residual']),abs(v['number_daughter_residual'])) for v in hist),
            scope='Hartree quasiparticle production, optional elastic, no zero mode or chemical reactions; collision Jacobian frozen gap; occupation guard is heuristic')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=64);p.add_argument('--sigma',type=float,default=1e-7)
    p.add_argument('--elastic',action='store_true');p.add_argument('--qmin',type=float,default=1e-5)
    p.add_argument('--end',type=float,default=.5);p.add_argument('--guard',type=float,default=1.)
    p.add_argument('--rtol',type=float,default=1e-5);p.add_argument('--wall',type=float,default=300.)
    p.add_argument('--output',required=True);a=p.parse_args()
    m=HartreePortal(sigma=a.sigma,elastic=a.elastic,table=PortalTable('.runtime/mesh_imPiR.dat'),n=a.n,nparent=80,kappa=.00685,
        scalar_ratio=1e-12,fermion_ratio=.9,qmin=a.qmin,g_fermion=4,parent_grid='balanced',fermion_grid='balanced',parent_floor=1e-18)
    r=m.run(a.end,a.rtol,a.wall,a.guard)
    Path(a.output).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print('FINAL',json.dumps(r['final']),r['message'],flush=True)
