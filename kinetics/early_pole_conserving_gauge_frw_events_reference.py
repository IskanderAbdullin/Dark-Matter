"""Fixed hard-support gauge mean field in an independently audited FRW toy."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from early_pole_conserving_gauge_events_reference import Model


class Expanding(Model):
    def __init__(self,g):
        super().__init__(g)
        self.node=np.concatenate(self.p)
        self.m=np.concatenate([np.full_like(p,m) for p,m in zip(self.p,self.M,strict=True)])
        self.rho_map=np.array([[0,1,1],[1,1,0],[1,0,1]])*(g*g/np.pi**4)

    def quantities(self,eta,z):
        a=np.exp(eta)
        a3=np.exp(-3*eta)
        p=self.node/a
        E=np.hypot(p,self.m)
        integrals=np.zeros((3,len(p)))
        for j,k in enumerate(self.labels[:3]):
            q=self.slices[k]
            integrals[j,q]=2*np.pi**2/(self.deg[j]*E[q])
        J=a3*integrals@z
        e=self.rho_map@J
        eg=e@integrals
        hf=.5*J@e
        r0=a3*(E@z)
        R0=a3*((self.m*self.m/E)@z)
        RK=-a3*(eg*self.m*self.m/(E*E))@z
        rho=r0+hf
        trace=R0+RK
        return dict(energy=E,gradient=E+eg,rho=rho,trace=trace,trace_field=RK,field=hf)

    def rates_time(self,eta,n,balanced):
        box=self.quantities(eta,n)
        f=n/self.w
        assert min(f)>0. and min(1+self.stat*f)>0.
        f0,f1=np.log(f),np.log1p(self.stat*f)
        d=np.zeros_like(n)
        prod=0.
        for v,k in self.events(box['gradient'] if balanced else box['energy']):
            p,m=np.maximum(v,0.),np.maximum(-v,0.)
            a=m@f0+p@f1
            b=p@f0+m@f1
            gain=k*np.exp(b)*np.expm1(a-b)
            d+=v*gain
            prod+=gain*(a-b)
        return d,prod,box

    def evolution(self,balanced):
        A=self.init.copy()
        e0=self.quantities(0.,A)
        h=3e-5
        def diff(key):
            a,b,c,d=(self.quantities(j*h,A)[key] for j in (-2,-1,1,2))
            return (a-8*b+8*c-d)/(12*h)
        partial=diff('rho')
        field_partial=diff('field')
        resid_field=field_partial+4*e0['field']-e0['trace_field']
        assert abs(resid_field)<2e-9*max(e0['field'],1e-40)
        assert abs(partial+4*e0['rho']-e0['trace'])<2e-8*e0['rho']
        def rhs(eta,w):
            d,prod,b=self.rates_time(eta,w[:-2],balanced)
            coll=np.exp(-3*eta)*b['gradient']@d
            return np.r_[d,np.exp(4*eta)*b['trace'],np.exp(4*eta)*coll]
        sol=solve_ivp(rhs,(0.,.025),np.r_[A,0.,0.],method='DOP853',rtol=3e-12,atol=1e-14)
        assert sol.success
        eta=sol.t[-1]; q=sol.y[:-2,-1]
        b=self.quantities(eta,q)
        trace_integral,coll=sol.y[-2:,-1]
        residual=np.exp(4*eta)*b['rho']-e0['rho']-trace_integral
        assert abs(residual-coll)<4e-13
        assert abs(self.charge@(q-A))<1e-13
        if balanced:
            assert abs(residual)<1e-13
        return dict(balanced=balanced,degenerate_particle_count_relative_change=self.charge@(q-A),
                    continuity_residual=residual,collision_work_integral=coll,
                    initial_visible_independent_rho=e0['rho'],
                    initial_field_energy=e0['field'],initial_field_trace=e0['trace_field'],
                    pressure_check_at_fixed_comoving_state=partial+4*e0['rho']-e0['trace'],
                    pressure_check_field_only=resid_field)


def run():
    return dict(scope='Finite, already-hard tagged comoving supports; fixed physical bare masses; '
                      'arbitrary event strength per log scale factor; no Hubble fit, exact momentum '
                      'angle, late spectrum or quantum medium rate.',
                rows=[dict(g=g,cases=[Expanding(g).evolution(b) for b in (False,True)])
                      for g in (0.,.005,.05)])


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
