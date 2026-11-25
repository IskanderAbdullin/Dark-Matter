"""Fixed thermal mother, massless daughters; direct daughter-energy FV kernel.
All energies in T, time in inverse rest width. Not an FRW prediction.
"""
import argparse
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import expit, roots_legendre
from pathlib import Path


class DaughterEnergyFV:
    def __init__(self, n=200, mass=.1, qmin=1e-5, qmax=32., order=12,
                 bath='integrated'):
        self.n = n
        self.mass = mass
        self.edges = np.geomspace(qmin, qmax, n+1)
        l, h = self.edges[:-1], self.edges[1:]
        self.vol = (h**3-l**3)/(6*np.pi**2)
        self.mean = .75*(h**4-l**4)/(h**3-l**3)
        self.mean2 = .6*(h**5-l**5)/(h**3-l**3)
        c = mass**2/4
        # q*r>=c, with q=E_phi/T, r=E_chi/T. Split at the two
        # intersections with the horizontal cell edges.
        x, w = roots_legendre(order)
        pref = 4*mass/(2*np.pi**2)
        qlo, qhi = l[:, None], h[:, None]
        rlo, rhi = l[None, :], h[None, :]
        start = np.maximum(qlo, c/rhi)
        split = np.clip(c/rlo, qlo, qhi)
        R = np.zeros((n,n)); A = np.zeros_like(R)
        for lo, hi in [(start, split), (np.maximum(qlo, c/rlo), qhi)]:
            width = np.maximum(hi-lo, 0.)
            for xi, wi in zip(x,w):
                q = lo + width*(xi+1)/2
                lower = np.minimum(rhi, np.maximum(rlo,c/q))
                dr = rhi-lower
                weight = wi*width/2
                R += weight*dr
                # Integral of Fermi function in r, using a stable log ratio.
                A += weight*np.log1p(np.exp(-q-lower)*(-np.expm1(-dr))
                                     /(1+np.exp(-q-rhi)))
        # Exact symmetry of the two energy coordinates; remove quadrature asymmetry.
        self.R = pref*(R+R.T)/2
        self.A = pref*(A+A.T)/2
        if bath == 'balanced':
            self.A = self.R*expit(-self.mean[:,None]-self.mean[None,:])
        self.asum = self.A.sum(axis=1)
        self.bath = bath

    def rhs(self, t, y):
        p, f = np.split(y, 2)
        return np.r_[(self.A@(1-f)+p*(self.asum-self.R@f))/self.vol,
                     (self.A@(1+p)-f*(self.asum+self.R@p))/(4*self.vol)]

    def jac(self, t, y):
        p, f = np.split(y, 2)
        v = self.vol
        return np.block([[np.diag((self.asum-self.R@f)/v),
                          -(self.A+p[:,None]*self.R)/v[:,None]],
                         [(self.A-f[:,None]*self.R)/(4*v[:,None]),
                          -np.diag((self.asum+self.R@p)/(4*v))]])

    def snapshot(self, y):
        p, f = np.split(y,2); weights = p*self.vol; N = weights.sum()
        # Exact fraction of phase volume below threshold for piecewise-constant f.
        l,h = self.edges[:-1], self.edges[1:]
        cold = np.maximum(np.minimum(h,.1)**3-l**3,0)/(h**3-l**3)
        return dict(number=float(N), mean_q=float(weights@self.mean/N),
                    rms_q=float(np.sqrt(weights@self.mean2/N)),
                    hot_fraction_q_gt_0p1=float(weights@(1-cold)/N),
                    peak_bin_q=float(self.mean[np.argmax(p)]),
                    max_occupation=float(max(p)), max_fermion=float(max(f)),
                    min_occupation=float(min(min(p),min(f))),
                    number_mismatch=float((N-(4*f*self.vol).sum())/N))

    def run(self, tend=.03, rtol=2e-7, max_step=.0003):
        # Independent vacuum gain integral, with q,r truncated by domain.
        lo,hi = self.edges[0],self.edges[-1]; c=self.mass**2/4
        def gain(q):
            low=max(lo,c/q)
            if low>=hi: return 0.
            return np.log1p(np.exp(-q-low)*(-np.expm1(-(hi-low)))/(1+np.exp(-q-hi)))
        exact=4*self.mass/(2*np.pi**2)*quad(gain,max(lo,c/hi),hi,
                              epsabs=1e-11,epsrel=1e-10,limit=200)[0]
        err=float(abs(self.A.sum()/exact-1))
        assert np.min(self.A)>=0 and np.max(self.A-self.R)<=1e-15
        if self.bath=='integrated': assert err<1e-6, err
        eq=np.r_[1/np.expm1(self.mean),expit(-self.mean)]
        balance=float(np.sum(np.abs(self.rhs(0,eq)[:self.n])*self.vol)/exact)
        if self.bath=='balanced': assert balance<1e-12, balance
        sol=solve_ivp(self.rhs,[0,tend],np.zeros(2*self.n),method='BDF',
                      jac=self.jac,t_eval=[tend],rtol=rtol,atol=1e-12,max_step=max_step)
        if not sol.success: raise RuntimeError(sol.message)
        snap=self.snapshot(sol.y[:,-1]); snap.update(n=self.n,bath=self.bath,
            qmin=lo,qmax=hi,tend=tend,rtol=rtol,max_step=max_step,
            spontaneous_gain_relative_error=err,discrete_thermal_residual=balance,nfev=sol.nfev)
        assert abs(snap['number_mismatch'])<1e-8 and snap['min_occupation']>-1e-8
        assert snap['max_fermion']<1.00001
        return snap


def main():
    p=argparse.ArgumentParser(); p.add_argument('--n',nargs='+',type=int,default=[100,200,400]);
    p.add_argument('--bath',choices=['integrated','balanced'],default='integrated');
    p.add_argument('--qmin',type=float,default=1e-5);p.add_argument('--qmax',type=float,default=32.)
    p.add_argument('--rtol',type=float,default=2e-7);p.add_argument('--max-step',type=float,default=.0003)
    p.add_argument('--tend',type=float,default=.03);p.add_argument('--output',required=True)
    a=p.parse_args(); rows=[]
    for n in a.n:
        result=DaughterEnergyFV(n,qmin=a.qmin,qmax=a.qmax,bath=a.bath).run(a.tend,a.rtol,a.max_step)
        rows.append(result);print(json.dumps(result),flush=True)
    Path(a.output).write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__': main()
