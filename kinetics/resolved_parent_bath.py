"""Momentum-resolved mother and daughters, with conservative dark events.

Vacuum isotropic 1<->2 decay; constant SM g; phenomenological bath rates.
The two-daughter allowed FV measure is projected onto adjacent parent
representative energies. Fractional stoichiometry conserves represented
particle number and energy. Geometric interpolation of both occupancies
and blocking factors gives detailed balance and inward boundary fluxes.
A partial cell whose representative pair energy is outside the parent
grid is discarded and its equilibrium-weighted measure is reported.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre,expit
from scipy.integrate import BDF,Radau
from scipy.optimize import brentq
from types import SimpleNamespace
from two_mass_bath import TwoMassBath
from radiation_fixed_bath import omega_from_number
from yukawa_rate_bath import rate_mapping


class ResolvedParent(TwoMassBath):
    def __init__(self,n=120,nparent=120,kappa_sm1=47.,kappa2=.03,gamma=.01,
                 scatter=True,parent_floor=0.,**kwargs):
        super().__init__(n=n,**kwargs)
        assert nparent>=32 and self.fermion_ratio>0 and kappa_sm1>=0 and kappa2>=0
        self.n1=nparent;self.dim=2*n+nparent
        self.kappa_sm1=kappa_sm1;self.kappa2=kappa2;self.gamma=gamma;self.scatter=scatter
        assert 0<=parent_floor<=1e-8
        self.parent_floor=parent_floor
        # Parent momentum support covers both daughter momenta, up to 2*qmax.
        nl=max(4,nparent//8);nm=3*nparent//4;nh=nparent-nl-nm
        qmax=max(90.,2*self.edges[-1]+10.)
        self.parent_edges=np.r_[0.,np.geomspace(1e-4,.2,nl)[0:],
            np.linspace(.2,24.,nm+1)[1:],np.geomspace(24.,qmax,nh+1)[1:]]
        assert len(self.parent_edges)==nparent+1
        l,h=self.parent_edges[:-1],self.parent_edges[1:]
        self.v1=self.g_fermion*(h-l)*(h*h+h*l+l*l)/(6*np.pi**2)
        u,w=roots_legendre(8)
        self.p1=(h[:,None]+l[:,None])/2+(h-l)[:,None]*u/2
        self.w1=(h-l)[:,None]/2*w*self.p1**2
        self.w1/=self.w1.sum(axis=1)[:,None]
        self.v2=self.g_fermion*self.vol
        self.vols=np.r_[self.vol,self.v2,self.v1]
        self.source1=np.r_[np.zeros(2*n),self.v1]
        self.source2=np.r_[np.zeros(n),self.v2,np.zeros(nparent)]
        self.max_dropped_equilibrium_fraction=0.
        self.max_projection_energy_error=0.

    def geometry(self,x):
        if x==self.geom_x:return
        super().geometry(x)
        e=np.hypot(self.p1,x)
        self.e1=np.sum(self.w1*e,axis=1)
        self.t1=np.sum(self.w1*x*x/e,axis=1)
        self.energies=np.r_[self.energy_phi,self.energy_chi,self.e1]
        self.energy_weights=self.vols*self.energies
        self.traces=np.r_[self.trace_weights,self.v1*self.t1]
        i,j=np.nonzero(self.R0>0)
        E=self.pair_E[i,j];R=self.R0[i,j]
        valid=(E>=self.e1[0])&(E<=self.e1[-1])
        eqweighted=R*expit(-E)
        loss=float(eqweighted[~valid].sum()/max(eqweighted.sum(),1e-300))
        self.max_dropped_equilibrium_fraction=max(self.max_dropped_equilibrium_fraction,loss)
        i,j,E,R=i[valid],j[valid],E[valid],R[valid]
        hi=np.clip(np.searchsorted(self.e1,E,side='right'),1,self.n1-1)
        lo=hi-1
        wh=np.clip((E-self.e1[lo])/(self.e1[hi]-self.e1[lo]),0.,1.)
        wl=1-wh
        self.max_projection_energy_error=max(self.max_projection_energy_error,
            float(np.max(np.abs(wl*self.e1[lo]+wh*self.e1[hi]-E)/E,initial=0.)))
        self.ii,self.jj,self.lo,self.hi=i,j,lo,hi
        self.wl,self.wh,self.event_rate0=wl,wh,R
        self.event_idx=np.vstack((i,self.n+j,2*self.n+lo,2*self.n+hi))
        self.stoich=np.vstack((1/self.vol[i],1/self.v2[j],-wl/self.v1[lo],-wh/self.v1[hi]))
        self.jidx=(self.event_idx[:,None,:]*self.dim+self.event_idx[None,:,:]).ravel()

    def geometric(self,a,b,wl,wh,jac=False):
        # Exact weighted powers at positive occupations. Optional C1
        # continuation below a numerical floor removes the fractional-power
        # derivative singularity without changing ordinary thermal bins.
        # At zero the forward availability remains exactly zero.
        def power(v,w):
            floor=self.parent_floor
            vv=np.maximum(v,1e-280)
            val=np.exp(w*np.log(vv))
            der=val*w/vv if jac else None
            if floor>0:
                low=(v<floor)&(w>0)
                t=np.clip(v[low]/floor,0.,1.)
                u=w[low];scale=np.exp(u*np.log(floor))
                val[low]=scale*((2-u)*t-(1-u)*t*t)
                if jac:
                    der[low]=np.where(v[low]>=0,scale*((2-u)-2*(1-u)*t)/floor,0.)
            else:
                bad=(v<=0)&(w>0)
                val[bad]=0.
                if jac:der[bad]=0.
            val[w==0]=1.
            if jac:der[w==0]=0.
            return val,der
        A,dA=power(a,wl);B,dB=power(b,wh)
        if jac:return A*B,dA*B,A*dB
        return A*B

    def dark(self,yd,Hf=1.,jac=False):
        n=self.n
        phi=yd[self.ii];f=yd[n+self.jj]
        l=yd[2*n+self.lo];h=yd[2*n+self.hi]
        # Parent boundary flux remains inward even in fractional events.
        if jac:
            F,Fl,Fh=self.geometric(l,h,self.wl,self.wh,True)
            B,Bl,Bh=self.geometric(1-l,1-h,self.wl,self.wh,True)
        else:
            F=self.geometric(l,h,self.wl,self.wh)
            B=self.geometric(1-l,1-h,self.wl,self.wh)
        R=self.event_rate0*Hf
        forward=F*(1-f);inverse=B*f
        net=R*(forward+phi*(forward-inverse))
        dy=np.bincount(self.event_idx.ravel(),weights=(self.stoich*net).ravel(),minlength=self.dim)
        if not jac:return dy
        deriv=np.vstack((R*(forward-inverse),-R*(F*(1+phi)+B*phi),
            R*((1+phi)*(1-f)*Fl+phi*f*Bl),
            R*((1+phi)*(1-f)*Fh+phi*f*Bh)))
        J=np.bincount(self.jidx,weights=(self.stoich[:,None,:]*deriv[None,:,:]).ravel(),
                      minlength=self.dim*self.dim).reshape(self.dim,self.dim)
        return dy,J

    def rates(self,x,theta,Hf):
        b=x*self.fermion_ratio
        r1=self.kappa_sm1*x*x*Hf*self.t1/x
        r2=self.kappa2*x*x*Hf*self.trace_weights[self.n:]/(self.v2*b)
        t1=np.zeros_like(r1);t2=np.zeros_like(r2)
        if self.scatter:
            U=x/theta;V=b/theta
            sc1=8*np.pi*self.gamma*self.kappa_sm1*x*theta*Hf/(1+U*U)
            sc2=8*np.pi*self.gamma*self.kappa2*x*theta*Hf/(self.fermion_ratio*(1+V*V))
            r1+=sc1;r2+=sc2;t1+=sc1*(1+2*U*U/(1+U*U));t2+=sc2*(1+2*V*V/(1+V*V))
        return np.r_[np.zeros(self.n),r2,r1],np.r_[np.zeros(self.n),t2,t1]

    def evaluate(self,z,y,jac=False):
        x=np.exp(z);self.geometry(x)
        yd=y[:self.dim];theta=np.exp(y[-4]);rb=self.rad0*theta**4
        total=rb+self.energy_weights@yd;Hf=np.sqrt(self.rad0/total)
        rate,rateT=self.rates(x,theta,Hf)
        eq=expit(-self.energies/theta)
        bath=rate*(eq-yd)
        if jac:
            dark,J=self.dark(yd,Hf,True)
        else:dark=self.dark(yd,Hf)
        dy=dark+bath
        # The SM receives bath exchange ONLY; dark events must conserve energy.
        dt=-self.energy_weights@bath/(4*rb)
        rhs=np.r_[dy,dt,self.traces@yd,self.source1@bath,self.source2@bath]
        if not jac:return rhs
        dH=-.5/total
        J-=np.diag(rate)
        J+=dH*np.outer(dy,self.energy_weights)
        bathT=rate*eq*(1-eq)*self.energies/theta+rateT*(eq-yd)
        JT=bathT+dH*dy*4*rb
        B=-np.diag(rate)+dH*np.outer(bath,self.energy_weights)
        BT=bathT+dH*bath*4*rb
        jtemp=-self.energy_weights@B/(4*rb)
        tempT=-self.energy_weights@BT/(4*rb)-4*dt
        D=np.column_stack((J,JT,np.zeros((self.dim,3))))
        return np.vstack((D,np.r_[jtemp,tempT,0.,0.,0.],np.r_[self.traces,0.,0.,0.,0.],
            np.r_[self.source1@B,self.source1@BT,0.,0.,0.],
            np.r_[self.source2@B,self.source2@BT,0.,0.,0.]))

    def rhs(self,z,y):return self.evaluate(z,y)
    def jac(self,z,y):return self.evaluate(z,y,True)

    def validate(self):
        checks={}
        self.geometry(.2);eq=expit(-self.energies);eq[:self.n]=1/np.expm1(self.energy_phi)
        c=self.dark(eq)
        balance=float(np.abs(c)@self.vols/max(self.event_rate0@expit(-self.pair_E[self.ii,self.jj]),1e-300))
        assert balance<1e-8,balance
        rng=np.random.default_rng(47)
        y=np.r_[rng.uniform(.05,.3,self.n),rng.uniform(.03,.2,self.n+self.n1),.01,0.,0.,0.]
        c=self.dark(y[:self.dim])
        scale=np.abs(c)@self.energy_weights
        energy=float(abs(c@self.energy_weights)/max(scale,1e-300))
        numbers=[float(self.vol@c[:self.n]),float(self.v2@c[self.n:2*self.n]),float(self.v1@c[2*self.n:])]
        num=max(abs(numbers[0]-numbers[1]),abs(numbers[0]+numbers[2]))/max(np.abs(c)@self.vols,1e-300)
        assert energy<1e-12 and num<1e-12,(energy,num)
        dr=rng.normal(size=len(y));dr/=np.linalg.norm(dr);eps=1e-6
        fd=(self.rhs(np.log(.2),y+eps*dr)-self.rhs(np.log(.2),y-eps*dr))/(2*eps)
        jd=self.jac(np.log(.2),y)@dr
        jerr=float(np.linalg.norm(fd-jd)/max(np.linalg.norm(jd),1e-300));assert jerr<1e-5,jerr
        for bound in (0.,1.):
            q=y[:self.dim].copy();q[2*self.n::3]=bound;q[self.n:2*self.n:3]=bound
            dc=self.dark(q)
            mask=np.r_[np.zeros(self.n,dtype=bool),q[self.n:]==bound]
            assert np.min(dc[mask]*(1-2*bound),initial=0.)>=-1e-10
        gain=[]
        for x in [.1,.2,1.,5.,20.]:
            self.geometry(x);p=expit(-self.e1)
            c=self.dark(np.r_[np.zeros(2*self.n),p])
            expected=self.kappa*x*x*(self.v1*self.t1/x)@p
            actual=self.vol@c[:self.n]
            gain.append(dict(x=x,ratio=float(actual/expected)))
        checks.update(balance=balance,energy=energy,number=float(num),jacobian=jerr,
            spontaneous_rate=gain,max_projection_energy_error=self.max_projection_energy_error,
            max_dropped_equilibrium_fraction=self.max_dropped_equilibrium_fraction)
        assert max(abs(r['ratio']-1) for r in gain)<.15,gain
        return checks

    def evolve(self,x0=.1,xend=80.,rtol=2e-6,initial_parent='empty',keep_psd=False,parent_gev=86000.,wall_limit=0.,fermion_atol=1e-15,method="BDF"):
        check=self.validate();print('checks',json.dumps(check),flush=True)
        self.geometry(x0);ini=expit(-self.e1) if initial_parent=='thermal' else np.zeros(self.n1)
        y0=np.r_[np.zeros(2*self.n),ini,0.,0.,0.,0.]
        init=float(self.v1@ini);rho0=float(self.rad0+(self.v1*self.e1)@ini)
        ts=sorted(set([x0]+[x for x in [.02,.03,.05,.075,.1,.125,.15,.175,.2,.3,.5,1,2,3,5,10,15,20,30,40,60,80,100,120] if x0<x<xend]+[xend]))
        atol=np.r_[np.full(self.n,1e-12),np.full(self.n+self.n1,fermion_atol),1e-13,1e-12,1e-13,1e-13]
        start=time.monotonic();last_print=start
        zlist=np.log(ts);ys=[y0];actual_z=[zlist[0]];k=1;steps=0;message='Integration complete.'
        solver={"BDF":BDF,"Radau":Radau}[method](self.rhs,zlist[0],y0,zlist[-1],jac=self.jac,rtol=rtol,atol=atol,max_step=.035)
        while solver.status=='running':
            err=solver.step();steps+=1
            if solver.status=='failed':
                message=str(err);break
            now=time.monotonic()
            if now-last_print>=15. or (k<len(zlist) and solver.t>=zlist[k]):
                print('progress',json.dumps(dict(x=float(np.exp(solver.t)),steps=steps,nfev=solver.nfev,
                    elapsed=now-start)),flush=True);last_print=now
            if k<len(zlist) and solver.t>=zlist[k]:
                interp=solver.dense_output()
                while k<len(zlist) and solver.t>=zlist[k]:
                    ys.append(interp(zlist[k]));actual_z.append(zlist[k]);k+=1
            if wall_limit>0 and now-start>=wall_limit and solver.status=='running':
                message='Explicit wall-clock budget reached.';break
        success=solver.status=='finished'
        if not success and solver.t>actual_z[-1]:
            actual_z.append(solver.t);ys.append(solver.y.copy())
        sol=SimpleNamespace(t=np.array(actual_z),y=np.array(ys).T,success=success,message=message,
            nfev=solver.nfev,njev=solver.njev)
        history=[]
        for x,y in zip(np.exp(sol.t[1:]),sol.y.T[1:]):
            self.geometry(x);p=y[:self.n];f=y[self.n:2*self.n];g=y[2*self.n:self.dim]
            zt,W,S1,S2=y[-4:];theta=np.exp(zt);s=self.snapshot(y[:2*self.n])
            Np=s['number'];N1=float(self.v1@g);N2=float(self.v2@f)
            r1=float((self.v1*self.e1)@g);r2=float((self.v2*self.energy_chi)@f)
            rp=float((self.vol*self.energy_phi)@p);total=self.rad0*theta**4+r1+r2+rp
            err1=(N1+Np-init-S1)/max(Np,1e-12);err2=(N2-Np-S2)/max(Np,1e-12)
            cont=(total-rho0-W)/rho0
            l,h=self.edges[:-1],self.edges[1:];cold=np.maximum(np.minimum(h,.1*theta)**3-l**3,0)/(h**3-l**3)
            s.pop('number_mismatch')
            s.update(x=x,theta=float(theta),mother_number=N1,mother_eq_number=float(self.v1@expit(-self.e1/theta)),
                mother_rho_Tref4=r1,fermion_number=N2,rho_chi2_Tref4=r2,rho_phi_Tref4=rp,total_rho_Tref4=float(total),
                mother_max_occupation=float(g.max()),mother_min_occupation=float(g.min()),integrated_parent_bath_source=float(S1),
                integrated_fermion_bath_source=float(S2),number_source_parent_residual=float(err1),number_source_residual=float(err2),
                continuity_residual=float(cont),matter_fraction_chi2=float(r2/total),
                mean_p_over_actual_T=float(s['mean_q']/theta),hot_fraction_p_over_actual_T_gt_0p1=float((p*self.vol)@(1-cold)/Np),
                omega_phi_h2_m86=float(omega_from_number(Np/theta**3,gs=self.g_bath)))
            assert abs(cont)<5e-5 and max(abs(err1),abs(err2))<1e-5,s
            assert s['min_occupation']>-1e-7 and s['max_fermion']<1.00001 and g.min()>-1e-7 and g.max()<1.00001,s
            if N1>1e-14:
                xi=brentq(lambda u:self.v1@expit(u-self.e1/theta)-N1,-750.,100.)
                fit=expit(xi-self.e1/theta)
                s['parent_shape_TV']=float(.5*self.v1@np.abs(g-fit)/N1)
                s['parent_energy_rel_difference']=float(r1/((self.v1*self.e1)@fit)-1)
                d=self.dark(y[:self.dim]);dc=self.dark(np.r_[p,f,fit])
                scale=max(abs(self.vol@d[:self.n]),abs(self.vol@dc[:self.n]),1e-30)
                s['closure_dark_source_difference']=float((self.vol@(dc[:self.n]-d[:self.n]))/scale)
            history.append(s)
            print('state',json.dumps({k:s[k] for k in ('x','omega_phi_h2_m86','hot_fraction_p_over_actual_T_gt_0p1','mother_number','mother_min_occupation','continuity_residual')}),flush=True)
        out=dict(n=self.n,nparent=self.n1,x0=x0,xend=float(np.exp(sol.t[-1])),requested_xend=xend,kappa=self.kappa,kappa2=self.kappa2,
            kappa_sm1=self.kappa_sm1,gamma=self.gamma,scatter=self.scatter,scalar_ratio=self.scalar_ratio,fermion_ratio=self.fermion_ratio,
            parent_GeV=parent_gev,g_bath=self.g_bath,g_fermion=self.g_fermion,rtol=rtol,initial_parent=initial_parent,initial_mother_number=init,
            initial_total_rho_Tref4=rho0,checks=check,method=method,fermion_atol=fermion_atol,nfev=sol.nfev,njev=sol.njev,steps=steps,wall_seconds=time.monotonic()-start,parent_floor=self.parent_floor,success=bool(sol.success),message=sol.message,
            max_dropped_equilibrium_fraction=self.max_dropped_equilibrium_fraction,max_projection_energy_error=self.max_projection_energy_error,
            history=history,final=history[-1] if history else None)
        if keep_psd:
            out['frozen_cell_spectrum']=dict(edges=self.edges.tolist(),f_phi=sol.y[:self.n,-1].tolist(),f_chi2=sol.y[self.n:2*self.n,-1].tolist(),
                parent_edges=self.parent_edges.tolist(),f_parent=sol.y[2*self.n:self.dim,-1].tolist(),g_phi=1,g_chi2=self.g_fermion,g_parent=self.g_fermion,
                theta_freezeout=history[-1]['theta'] if history else 1.)
        if sol.success:out['xend']=xend
        if self.kappa_sm1>0 and self.kappa2>0:
            out['mapping']=rate_mapping(parent_gev,self.fermion_ratio,self.kappa,self.kappa2,self.kappa_sm1)
        return out


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=120);p.add_argument('--nparent',type=int,default=120);p.add_argument('--g-fermion',type=int,choices=[2,4],default=4)
    p.add_argument('--kappa',type=float,default=.005278180162600994);p.add_argument('--kappa-sm1',type=float,default=47.)
    p.add_argument('--kappa2',type=float,default=.03);p.add_argument('--gamma',type=float,default=.01)
    p.add_argument('--fermion-ratio',type=float,default=.9);p.add_argument('--parent-gev',type=float,default=86000.)
    p.add_argument('--x0',type=float,default=.1);p.add_argument('--xend',type=float,default=80.)
    p.add_argument('--qmin',type=float,default=1e-8);p.add_argument('--rtol',type=float,default=2e-6)
    p.add_argument('--initial-parent',choices=['empty','thermal'],default='empty');p.add_argument('--no-scatter',action='store_true')
    p.add_argument('--parent-floor',type=float,default=0.);p.add_argument('--wall-limit',type=float,default=0.)
    p.add_argument('--fermion-atol',type=float,default=1e-15);p.add_argument('--method',choices=['BDF','Radau'],default='BDF')
    p.add_argument('--test-only',action='store_true');p.add_argument('--keep-psd',action='store_true');p.add_argument('--output')
    a=p.parse_args()
    m=ResolvedParent(n=a.n,nparent=a.nparent,g_fermion=a.g_fermion,kappa=a.kappa,kappa2=a.kappa2,kappa_sm1=a.kappa_sm1,
        gamma=a.gamma,scatter=not a.no_scatter,parent_floor=a.parent_floor,scalar_ratio=86./(a.parent_gev*1e9),fermion_ratio=a.fermion_ratio,qmin=a.qmin)
    if a.test_only:print(json.dumps(m.validate(),indent=2));return
    assert a.output
    r=m.evolve(a.x0,a.xend,a.rtol,a.initial_parent,a.keep_psd,a.parent_gev,a.wall_limit,a.fermion_atol,a.method)
    Path(a.output).write_text(json.dumps([r],indent=2,allow_nan=False)+'\n')
    print('final',json.dumps({k:r[k] for k in ('n','nparent','success','message','nfev','max_dropped_equilibrium_fraction','final')}),flush=True)
    if not r['success']:raise RuntimeError(r['message'])


if __name__=='__main__':main()
