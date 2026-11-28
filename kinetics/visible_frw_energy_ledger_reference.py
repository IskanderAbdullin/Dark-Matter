"""Visible FRW temperature and energy ledger with signed external heat."""
import argparse
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre

from sm_dark_entropy_exchange_reference import StandardModel

_QX, _QW = roots_legendre(12)


class ConservativeEOS:
    """Same interpolated entropy, pressure integrated sharply across each knot."""
    def __init__(self, eos):
        self.base = eos
        self.lo, self.hi = eos.lo, eos.hi
        self.edge = eos.gs.x
        self.anchor = eos.state(np.exp(self.lo))[1]
        incr = self.primitive(self.edge[:-1], self.edge[1:])
        self.cumulative = np.r_[0., np.cumsum(incr)]

    def primitive(self, lo, hi):
        from sm_dark_entropy_exchange_reference import C
        w = (hi-lo)/2
        x = np.asarray((hi+lo)/2)[...,None]+np.asarray(w)[...,None]*_QX
        return 4*C/3*w*((self.base.gs(x)*np.exp(4*x))@_QW)

    def state(self, T):
        from sm_dark_entropy_exchange_reference import C
        x = float(np.log(T))
        assert self.lo<=x<=self.hi
        j = min(max(0,int(np.searchsorted(self.edge,x,side="right")-1)),len(self.edge)-2)
        s = 4*C/3*T**3*self.base.gs(x)
        cv = 4*C/3*T**3*(3*self.base.gs(x)+self.base.slope(x))
        p = self.anchor+self.cumulative[j]+self.primitive(self.edge[j],x)
        return float(T*s-p),float(p),float(s),float(cv)


class VisibleBath:
    def __init__(self, eos, T):
        self.eos = eos
        self.a = 1.
        self.T = float(T)
        self.E0 = eos.state(T)[0]
        self.trace = 0.
        self.portal = 0.
        self.entropy = 0.
        self.entropy_first_order = 0.
        self.s0 = eos.state(T)[2]
        self.kick_residual = 0.

    def invert(self, index, target):
        """Invert a monotone EOS density or entropy, without extrapolation."""
        T = self.T
        eq = self.eos.state(T)[index]
        if target == eq:
            return T
        a, b = (self.eos.lo, np.log(T)) if target < eq else (np.log(T), self.eos.hi)
        return float(np.exp(brentq(lambda l: self.eos.state(np.exp(l))[index]/target-1.,
                                  a, b, xtol=2e-15, rtol=9e-16)))

    def advance(self, a):
        """Entropy-preserving redshift plus independent pressure trace."""
        assert a >= self.a
        if a == self.a:
            return
        base = self.a
        loga = np.log(a/base)
        s_old = self.eos.state(self.T)[2]
        E_before = self.energy
        trace = 0.
        after = self.invert(2, s_old*(base/a)**3)
        # Trace integrand inherits EOS interpolation knots. Split in log(a)
        # at each crossed knot, without using energy differences as tracer.
        knots = self.eos.edge[(self.eos.edge>np.log(after))&(self.eos.edge<np.log(self.T))]
        split = sorted([0.,loga,*[np.log(s_old/self.eos.state(np.exp(k))[2])/3 for k in knots]])
        for l,r in zip(split,split[1:]):
            for z,w in zip(_QX,_QW):
                aj=base*np.exp(l+(r-l)*(1+z)/2)
                tj=self.invert(2,s_old*(base/aj)**3)
                rho,p,*_=self.eos.state(tj)
                trace+=w*aj**4*(rho-3*p)*(r-l)/2
        self.T = after
        self.a = a
        assert abs(self.energy-E_before-trace)<5e-12*self.E0
        self.trace += trace

    def exchange(self, dark_a4_energy):
        """Withdraw visible energy for all daughter momenta (negative returns heat)."""
        if dark_a4_energy == 0.:
            return
        before,_,s,_ = self.eos.state(self.T)
        target = before-dark_a4_energy/self.a**4
        if target <= 0:
            raise ValueError("portal energy withdrawal exceeds visible bath")
        t_new = self.invert(0, target)
        rho,_,s_new,_ = self.eos.state(t_new)
        dE = (rho-before)*self.a**4
        self.portal += dark_a4_energy
        self.entropy += (s_new-s)*self.a**3
        self.entropy_first_order -= dark_a4_energy/(self.a*self.T)
        self.kick_residual = max(self.kick_residual,abs(dE+dark_a4_energy)/self.E0)
        assert self.kick_residual<3e-13
        self.T = t_new

    @property
    def energy(self):
        return self.a**4*self.eos.state(self.T)[0]

    def summary(self):
        e_res = (self.energy-self.E0-self.trace+self.portal)/self.E0
        assert abs(e_res)<4e-11
        s_res=(self.eos.state(self.T)[2]*self.a**3-self.s0-self.entropy)/self.s0
        assert abs(s_res)<3e-13
        return dict(T_visible_GeV=self.T,
                    visible_pressure_trace_a4_energy=self.trace,
                    visible_portal_withdrawal_a4_energy=self.portal,
                    visible_comoving_entropy_change=self.entropy,
                    visible_comoving_entropy_change_first_order=self.entropy_first_order,
                    visible_entropy_ledger_defect_in_initial_entropy_units=s_res,
                    initial_visible_energy=self.E0,
                    final_visible_a4_energy=self.energy,
                    visible_energy_ledger_defect_in_initial_visible_units=e_res,
                    max_visible_kick_defect_in_initial_visible_units=self.kick_residual)


def test(eos,T=.1):
    """Signed manufactured exchange, checked against independent differential IVP."""
    extent=.10
    k0=eos.state(T)[0]
    def flow(x):
        return k0*(.12*np.cos(70*x)-.06)
    def rhs(x,y):
        t=y[0]
        rho,p,s,cv=eos.state(t)
        return [(-3*t*s-flow(x)*np.exp(-4*x))/cv,flow(x),np.exp(4*x)*(rho-3*p)]
    ref=solve_ivp(rhs,(0.,extent),[T,0.,0.],rtol=5e-13,atol=2e-17,max_step=.002)
    assert ref.success
    out=[]
    for n in (6,12,24):
        st=VisibleBath(eos,T)
        dx=extent/n
        for i in range(n):
            x=(i+.5)*dx
            st.advance(np.exp(x))
            st.exchange(dx*flow(x))
            st.advance(np.exp((i+1)*dx))
        r=st.summary()
        out.append(dict(steps=n,relative_final_T_defect=abs(st.T/ref.y[0,-1]-1),**r))
    assert out[-1]["relative_final_T_defect"]<out[0]["relative_final_T_defect"]/12
    s=VisibleBath(eos,T)
    s.advance(np.exp(extent))
    assert abs(eos.state(s.T)[2]*np.exp(3*extent)/eos.state(T)[2]-1)<3e-14
    return dict(manufactured_signed=out,reference_final_T_GeV=ref.y[0,-1],no_portal=s.summary())


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("eos")
    a=p.parse_args()
    print(json.dumps(test(ConservativeEOS(StandardModel(a.eos))),indent=2,allow_nan=False))
