"""Restore the frozen exterior self-log in a conditional Friedmann ledger."""
import gzip
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_legendre

from finite_mass_qed_live_log_exterior_reference import WithExterior, clip_cells
from finite_mass_qed_shell_log_reference import integrated_kernel


class Complete(WithExterior):
    def __init__(self,data,eta_end):
        super().__init__(data,.25)
        self.dark_exterior=[]
        self.low_free=[]
        x,w=roots_legendre(16)
        self.eta_end=eta_end
        for name,sl in zip(('parent','fermion'),self.r[1:],strict=True):
            e=np.asarray(data[f'endpoint_{name}_momentum_edges'])
            f=np.asarray(data[f'endpoint_{name}_occupations'])
            left,right=self.cells[sl.start].lo,self.cells[sl.stop-1].hi
            m=self.cells[sl.start].m
            edges=np.unique(np.r_[e[e>.25],.25,left,right])
            mid=(edges[1:]+edges[:-1])/2
            vals=f[np.searchsorted(e,mid)-1]*((mid<left)|(mid>right))
            self.dark_exterior.append((edges,vals,m))
            # The previous log stress reference cuts the hard functional exactly
            # at .25. Restore bare charged free energy below that scale only.
            pp=[];ww=[]
            for l,u,f0 in clip_cells(e,f,0.,.25):
                p=(l+u)/2+(u-l)*x/2
                pp.extend(p)
                ww.extend(f0*(u-l)/2*w*2*p*p/np.pi**2)
            self.low_free.append((np.asarray(pp),np.asarray(ww),m))
        flo,fhi=self.edge_sets[0][0],self.edge_sets[0][-1]
        def integrate(power):
            def f(p):
                return p**power/np.expm1(p+.3)
            first=quad(f,0.,flo,epsabs=1e-13)[0]
            last=quad(f,fhi,60.,epsabs=1e-13)[0]
            return first+last
        self.external_photon_J=integrate(1)
        self.external_photon_R=integrate(3)/np.pi**2
        self.add_self=False
        self.slog=None
        self.audit=None

    def slow_exterior_log(self,eta,n=12):
        E=T=0.
        for edges,U,m in self.dark_exterior:
            K,D=integrated_kernel(edges,np.exp(eta)*m,n=n)
            E+=U@K@U
            T+=U@D@U
        return -self.c*E/4,-self.c*T/4

    def prepare(self):
        # The spectator's pure self energy only depends on eta. Use a
        # single interpolating energy and its derivative to retain the ledger.
        M=8
        z=np.cos(np.pi*(np.arange(M)+.5)/M)
        x=(z+1)*self.eta_end/2
        f=[self.slow_exterior_log(eta)[0] for eta in x]
        self.slog=Chebyshev.fit(x,f,M-1,[0.,self.eta_end])
        self.slog_prime=self.slog.deriv()
        errors=[]
        for eta in np.linspace(0.,self.eta_end,5):
            r,t=self.slow_exterior_log(eta)
            d1=abs(self.slog(eta)-r)
            d2=abs(self.slog_prime(eta)-t)
            if max(d1,d2)>2e-20:
                raise ValueError('Exterior self energy or scale work did not reproduce.')
            errors.append([float(d1),float(d2)])
        self.audit=np.max(errors,axis=0).tolist()
        lo=np.asarray(self.slow_exterior_log(self.eta_end/2,12))
        hi=np.asarray(self.slow_exterior_log(self.eta_end/2,20))
        self.mesh_check=float(max(abs(lo/hi-1)))
        if self.mesh_check>1e-7:
            raise ValueError('Pure exterior log quadrature did not resolve the cutoff shape.')

    def spectators(self,eta):
        J,D,R,T=super().spectators(eta)
        J[0]=self.external_photon_J
        R+=self.external_photon_R
        for p,w,m in self.low_free:
            x=m*np.exp(eta)
            e=np.hypot(p,x)
            R+=w@e
            T+=w@(x*x/e)
        if self.add_self:
            R+=self.slog(eta)
            T+=self.slog_prime(eta)
        return J,D,R,T


def run_clock(M,conditional):
    rhoV=np.pi**2*106.75/30
    end=M.eta_end
    ref_time=.04
    alpha=ref_time/end
    y0=np.r_[M.init.ravel(),M.Ns,0.,0.,0.,0.,0.]
    M.log_weight=1.
    M.cross_scale=1.
    M.add_self=True
    M.reject=[0.,0.]
    M.pair_scale=1.
    def energy(y,eta):
        M.eta=eta
        R,S=M.weights(y[:-6].reshape(-1,2),True)[-2:]
        _,r,s,_=M.scalar(y[-6])
        return R+r,S+s,M.trace
    start,s0,_=energy(y0,0.)
    effects=[]
    def rhs(eta,y):
        r,_,_=energy(y,eta)
        if r<0.:
            raise ValueError('Negative archived dark energy in the tested domain.')
        c=1/np.sqrt(1+r/rhoV)
        ext=M.slog(eta)
        d=-ext/(rhoV+r)
        rel=d/(1+np.sqrt(1+d))
        effects.append((c,rel))
        d=M.rhs_joint(y[:-1],True,1.)
        return np.r_[alpha*(c if conditional else 1.)*d,M.trace]
    sol=solve_ivp(rhs,(0.,end),y0,method='DOP853',rtol=8e-11,atol=2e-13,max_step=end/5)
    if not sol.success:
        raise ValueError('Conditional exterior evolution did not finish.')
    y1=sol.y[:,-1]
    r,s,_=energy(y1,end)
    absolute_drift=abs(r-start-y1[-1]-y1[-4])
    entropy=abs(s-s0-y1[-5])
    counts=abs(np.sum((y1-y0)[:-6].reshape(-1,2)[:,0]))
    if max(absolute_drift,counts)>2e-13 or entropy>3e-11:
        raise ValueError('Full spectator ledger did not close.')
    pure_bound=abs(M.slog(0.))*np.exp(2*end)/(2*rhoV)
    if max(x[1] for x in effects)>pure_bound:
        raise ValueError('Pure-spectator clock exceeded its diagnostic bound.')
    return dict(conditional_dark_H=conditional,scalar_gain=float(y1[-6]-M.Ns),
                dark_relative_energy_initial=float(start/rhoV),dark_relative_energy_final=float(r/rhoV),
                redshift_work=float(y1[-1]),energy_ledger_defect=float(absolute_drift),
                entropy_production=float(s-s0),entropy_ledger_defect=float(entropy),
                sampled_clock_factor_range=[float(min(x[0] for x in effects)),float(max(x[0] for x in effects))],
                sampled_pure_log_clock_increase_max=float(max(x[1] for x in effects)),
                analytical_pure_log_clock_bound=float(pure_bound))


def run():
    rows=[]
    for name in ('fine','third'):
        archive=Path(__file__).parent/f'results/early-pole-spectrum-{name}-m06-2026-11-26.json.gz'
        with gzip.open(archive) as f:d=json.load(f)
        M=Complete(d,.012)
        M.prepare()
        rows.append(dict(archive=name,external_self_log_initial=float(M.slog(0.)),
                         external_self_log_final=float(M.slog(.012)),self_log_tabulation_defects=M.audit,
                         self_log_quadrature_relative_defect=M.mesh_check,
                         imposed_clock=run_clock(M,False),conditional_clock=run_clock(M,True)))
    return dict(scope='Finite tagged archived-plus-fixed spectator control. Optional conditional Friedmann clock, arbitrary reaction calibration, hard-only log.',
                rows=rows)


if __name__=='__main__':
    print(json.dumps(run(),indent=2,allow_nan=False))
