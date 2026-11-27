"""Positive two-node parent collision stencil with exact projected energy.

An independent discretization of the held-expansion, live-parent control.
"""

import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit

from early_pole_decay_mode_clock_reference import G, KAPPA, X
from early_pole_resolved_parent_reference import ResolvedParent, S


class EnergyStencil(ResolvedParent):
    def __init__(self,ns,nf,nparent):
        super().__init__(ns,nf,nparent,0.,False)
        self.ed=.75*np.diff(self.d.fed**4)/np.diff(self.d.fed**3)

    def parent_energies(self,x):
        return np.sum(self.pqw*np.hypot(x,self.pq),axis=1)/self.pw

    def daughter_mass(self,x):
        return 0.

    def daughter_energies(self,x):
        return self.ed

    def channels(self,x,m):
        E=np.hypot(m,self.d.p)
        P=self.parent_energies(x)
        mass=self.daughter_mass(x)
        if x<=m+mass:
            z=np.zeros(0,dtype=int)
            f=np.zeros(0,dtype=float)
            return E,P,z,z,z,z,f,f,f,0.
        D=np.sqrt((x*x-(m+mass)**2)*(x*x-(m-mass)**2))/2
        A=(x*x+m*m-mass*mass)/2
        lo=(x*x*E*E+D*D)/(A*E+self.d.p*D)
        hi=lo+2*self.d.p*D/(m*m)
        edges=np.hypot(self.d.ef,mass)
        a=np.maximum(lo[:,None],E[:,None]+edges[:-1])
        b=np.minimum(hi[:,None],E[:,None]+edges[1:])
        rate=KAPPA*G*x*((x+mass)**2-m*m)/(E*self.d.p)
        w=np.maximum(b-a,0.)*(rate*self.d.w)[:,None]
        i,j=np.nonzero(w)
        ep=E[i]+self.daughter_energies(x)[j]
        J=np.searchsorted(P,ep,side='right')
        keep=(J>0)&(J<len(P))
        # The test for omitted support is FD weighted, not an unqualified
        # claim that all kinematically possible nonthermal tails are covered.
        W=w[i,j]
        expect=W*expit(-ep)
        miss=float(expect[~keep].sum()/max(expect.sum(),1e-300))
        i,j,ep,J,W=(q[keep] for q in (i,j,ep,J,W))
        delta=(ep-P[J-1])/(P[J]-P[J-1])
        return E,P,i,j+self.ns,J-1+self.iP,J+self.iP,1.-delta,delta,W,miss

    @staticmethod
    def parent_products(u,v,a,b,calc_jac):
        u,v=np.clip(u,0.,1.),np.clip(v,0.,1.)
        F=np.power(u,a)*np.power(v,b)
        B=np.power(1-u,a)*np.power(1-v,b)
        if not calc_jac:
            return F,B
        # Derivative at the interior states used for the physical BDF
        # trajectory. At an exactly clipped endpoint use a one-sided
        # finite zero convention for this numerical preconditioner only.
        def div(num,den):
            return np.divide(num,den,out=np.zeros_like(num),where=den>1e-100)
        return F,B,div(a*F,u),div(b*F,v),-div(a*B,1-u),-div(b*B,1-v)

    def at(self,x,m,y,jac=False):
        E,P,i,j,k,l,a,b,W,_=self.channels(x,m)
        fs,fd=y[i],y[j]
        K=self.parent_products(y[k],y[l],a,b,jac)
        F,B=K[:2]
        emit=(1+fs)*(1-fd)
        absorb=fs*fd
        R=W*(F*emit-B*absorb)
        pieces=((i,1.),(j,1.),(k,-a),(l,-b))
        v=np.zeros(self.n)
        for ix,f in pieces:
            v+=np.bincount(ix,weights=f*R,minlength=self.n)
        v/=self.w
        if not jac:
            return v
        dFdk,dFdl,dBdk,dBdl=K[2:]
        inputs=((i,W*(F*(1-fd)-B*fd)),
                (j,-W*(F*(1+fs)+B*fs)),
                (k,W*(dFdk*emit-dBdk*absorb)),
                (l,W*(dFdl*emit-dBdl*absorb)))
        J=np.zeros((self.n,self.n))
        for dim,c in pieces:
            for src,proj in inputs:
                np.add.at(J,(dim,src),c*proj)
        return J/self.w[:,None]

    def stencil_checks(self,x=.008,m=8e-5):
        E,P,i,j,k,l,a,b,W,miss=self.channels(x,m)
        Ed=self.daughter_energies(x)
        eps=np.max(abs((a*P[k-self.iP]+b*P[l-self.iP])/
                       (E[i]+Ed[j-self.ns])-1))
        pt=expit(-P)
        nd=expit(-Ed)
        ns=1/np.expm1(E)
        F,B=self.parent_products(pt[k-self.iP],pt[l-self.iP],a,b,False)
        fo=F*(1+ns[i])*(1-nd[j-self.ns])
        re=B*ns[i]*nd[j-self.ns]
        bal=np.max(abs(fo/re-1))
        y=np.r_[ns*.035,nd*.6,pt*.8]
        v=self.at(x,m,y)
        er=np.r_[E,Ed,P]*self.w
        norm=sum(abs(er*v))
        e_closure=abs(er@v)/norm
        assert max(eps,bal,e_closure)<3e-14
        # Probe sign at one parent-cell edge, not all parents simultaneously.
        qi=int(k[np.argmax(W)])
        for sign in [0,1]:
            test=y.copy()
            test[qi]=sign
            assert (1-2*sign)*self.at(x,m,test)[qi]>=-1e-15
        dr=np.cos(np.arange(self.n))*y
        dy=2e-5
        fd=(self.at(x,m,y+dy*dr)-self.at(x,m,y-dy*dr))/(2*dy)
        an=self.at(x,m,y,True)@dr
        jac_err=np.max(abs(fd-an))/np.max(abs(fd))
        assert jac_err<1e-7
        return dict(event_energy_error_relative=float(eps),detailed_balance_error_relative=float(bal),
                    bulk_energy_error_relative=float(e_closure),jacobian_error_relative=float(jac_err),
                    omitted_thermal_seed_fraction=miss)

    def run(self,xend):
        checks=self.stencil_checks()
        initial=np.r_[np.zeros(self.iP),expit(-self.parent_energies(X))]
        tol=np.r_[np.full(self.ns,1e-7),np.full(self.nf,8e-12),
                  np.maximum(initial[self.iP:]*5e-10,1e-19)]
        sol=solve_ivp(self.rhs,(0.,np.log(xend/X)),initial,method='BDF',
                      jac=self.jac,rtol=2e-8,atol=tol,max_step=.02,dense_output=True)
        assert sol.success,sol.message
        et=np.linspace(0,sol.t[-1],33)
        out=sol.sol(et)
        n_s=self.d.w@out[:self.ns]
        n_d=self.d.wf@out[self.ns:self.iP]
        pn=self.pw@(out[self.iP:]-initial[self.iP:,None])
        invariant=max(np.max(abs(n_s-n_d)),np.max(abs(n_s+pn)))
        assert invariant<3e-12
        assert min(out[:self.ns].min(),0.)>-6e-7 and out[self.ns:].min()>-1e-9
        assert out[self.ns:].max()<=1+1e-9
        closure=[]
        for x,y in zip(X*np.exp(et),out.T,strict=True):
            m=self.d.gap(x,y[:self.ns])
            E,P,*_,miss=self.channels(x,m)
            rate=self.at(x,m,y)
            ec=self.w*np.r_[E,self.daughter_energies(x),P]
            norm=sum(abs(ec*rate))
            closure.append((abs(ec@rate)/max(norm,1e-300),miss))
        audits=np.array(closure)
        assert audits[:,0].max()<3e-13
        return dict(cells=[self.ns,self.nf,self.np],xend=xend,scalar_yield=float(n_s[-1]/S),
                    maximum_absolute_sampled_number_invariant=float(invariant),
                    maximum_relative_collision_energy_residual=float(audits[:,0].max()),
                    max_omitted_thermal_seed_fraction=float(audits[:,1].max()),
                    parent_initial_occupation=initial[self.iP:].tolist(),
                    parent_final_occupation=out[self.iP:,-1].tolist(),
                    parent_shell_edges=self.pe.tolist(),static_checks=checks,solver_calls=sol.nfev)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scalar',type=int,default=35)
    p.add_argument('--fermion',type=int,default=40)
    p.add_argument('--parent',type=int,default=40)
    p.add_argument('--xend',type=float,default=.05)
    p.add_argument('--static-only',action='store_true')
    a=p.parse_args()
    d=EnergyStencil(a.scalar,a.fermion,a.parent)
    print(json.dumps(d.stencil_checks() if a.static_only else d.run(a.xend),
                     indent=2,allow_nan=False))
