"""Couple a tunable cold 4-to-2 threshold comparator to the early decay system."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

from early_pole_decay_mode_clock_reference import X
from early_pole_live_parent_energy_reference import state_moments
from early_pole_massive_daughter_energy_stencil_reference import MassiveFermion
from early_pole_resolved_parent_reference import S
from early_pole_scalar_number_change_stencil_reference import NumberChange


class Chemical(MassiveFermion):
    def __init__(self,ns,nf,nparent,ratio,strength):
        super().__init__(ns,nf,nparent,ratio,0.,True)
        self.strength=strength
        self.qlo=1e-6

    def collider(self,x,m):
        d=NumberChange(self.d.p,self.d.ed,m,self.qlo,self.strength)
        h=np.sqrt(np.pi**2*106.75/90)*(86./x)/2.435e18
        return d,h

    def at(self,x,m,y,jac=False):
        v=super().at(x,m,y,jac)
        if self.strength==0 or m<=2.5e-4:
            return v
        c,h=self.collider(x,m)
        ix=c.active
        F=y[ix]
        if not jac:
            v[ix]+=c.flux(F)/(c.vol*h)
            return v
        # Analytic fixed-m derivative; the outer Hartree derivative is
        # evaluated by the base implicit solver using the whole RHS.
        f=np.maximum(F,1e-200)
        lf=np.log(f)
        A=np.logaddexp(0.,-lf)
        chem=c.S@A
        left=c.A@lf+A[c.incoming].sum(axis=1)
        fw=c.W*np.exp(left+chem)
        bw=c.W*np.exp(left)
        R,_=c.events(F)
        pre=c.A.multiply(R[:,None]).multiply((1./f)[None,:])
        din=(c.A-c.S)*.5
        dout=(c.A+c.S)*.5
        derivative=pre+(dout.multiply(fw[:,None])-din.multiply(bw[:,None])).multiply(
           (-1./(f*(1+f)))[None,:])
        j=(c.S.T@derivative).toarray()/(c.vol[:,None]*h)
        v[np.ix_(ix,ix)]+=j
        return v

    def chem_number_rate(self,t,q):
        x=X*np.exp(t)
        E,n,p,m=state_moments(self,x,q)
        if self.strength==0 or m<=2.5e-4:
            return 0.
        c,h=self.collider(x,m)
        rx,_=c.events(q[c.active])
        hcorrect=np.sqrt(self.radiation/(self.radiation+n+p))
        return float(-2*sum(rx)*hcorrect/h)

    def from_snapshot(self,s,xend):
        x0=s['x_end']
        t0=np.log(x0/X)
        tf=np.log(xend/x0)
        y0=np.r_[s['endpoint_scalar_occupations'],s['endpoint_fermion_occupations'],s['endpoint_parent_occupations']]
        assert np.array_equal(self.d.ed,s['scalar_momentum_edges'])
        eps=np.r_[np.full(self.ns,1e-8),np.full(self.nf,2e-12),
                  np.maximum(y0[self.iP:]*2e-10,1e-19)]
        def stop(t,y):
            return self.d.gap(x0*np.exp(t),y[:self.ns])-4e-4
        stop.terminal=True
        stop.direction=-1.
        sol=solve_ivp(lambda t,y:self.rhs(t0+t,y),(0.,tf),y0,method='BDF',
                      jac=lambda t,y:self.jac(t0+t,y),
                      atol=eps,rtol=2e-8,max_step=.006,dense_output=True,events=stop)
        assert sol.success,(sol.message,sol.t[-1],sol.nfev)
        _,D0,P0,_=state_moments(self,x0,y0)
        y=sol.y[:,-1]; xf=x0*np.exp(sol.t[-1])
        _,D1,P1,m=state_moments(self,xf,y)
        values=[]
        for order in [2,4]:
            Z,W=roots_legendre(order)
            h=np.diff(sol.t)/2
            ti=(sol.t[:-1,None]+h[:,None]*(Z+1)).ravel()
            wi=(h[:,None]*W).ravel()
            vals=np.array([[self.trace(t0+t,sol.sol(t)),self.chem_number_rate(t0+t,sol.sol(t))] for t in ti])
            values.append(wi@vals)
        n0=self.d.w@y0[:self.ns]; n1=self.d.w@y[:self.ns]
        pn=self.pw@(y[self.iP:]-y0[self.iP:])
        fn=self.d.wf@(y[self.ns:self.iP]-y0[self.ns:self.iP])
        er=((D1-D0)+(P1-P0)-values[-1][0])/D1
        chem=((n1-n0)-fn-values[-1][1])/n0
        assert abs(er)<2e-6 and abs(chem)<8e-6 and abs(pn+fn)<2e-12,(er,chem,values)
        return dict(x0=x0,xend=xf,requested_xend=xend,model_threshold_strength=self.strength,
                    scalar_yield_initial=n0/S,scalar_yield_final=n1/S,
                    net_chemical_count_change_relative_to_initial_scalar=float(values[-1][1]/n0),
                    net_decay_count_change_relative_to_initial_scalar=float(fn/n0),
                    final_scalar_pole_over_T=m,relative_energy_closure=float(er),
                    relative_chemical_number_closure=float(chem),
                    paused_before_invalid_cold_limit=bool(sol.t_events[0].size),
                    solver_rhs_calls=sol.nfev)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input')
    p.add_argument('--strength',type=float,default=1e-4)
    p.add_argument('--xend',type=float,default=.052)
    a=p.parse_args()
    f=Path(a.input)
    raw=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
    s=json.loads(raw)
    shape={(67,112,122):(65,90,72),(97,156,176):(95,125,108)}
    d=Chemical(*shape[tuple(s['cells'])],s['mass_ratio'],a.strength)
    print(json.dumps(d.from_snapshot(s,a.xend),indent=2,allow_nan=False))
