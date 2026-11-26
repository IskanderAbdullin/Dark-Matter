"""Positive, two-moment visible injection on finite comoving supports.

Count all light stable fermionic daughters together (four degrees of freedom).
Rates are physical per unit time and volume. A three-channel ledger retains
direct light, prompt heavy-to-light, and prompt heavy-to-scalar separately.
The visible source is vacuum final state (no dark factors or inverse portal);
heavy parents are assumed to decay immediately and isotropically.
"""
import argparse
import json

import numpy as np
from empty_lepton_pair_spectrum_reference import Source
from three_species_decay_onsager_reference import M, m, mp, P, Es
from sm_dark_entropy_exchange_reference import StandardModel, MPL


def kinetic(p, mass):
    return p*p/(np.hypot(p, mass)+mass)


def direct_intervals(T, source, edges, batch=2048):
    """Analytic angular integration of the direct source's l=0,2 weights."""
    q,p,v,N,A = source.event_moments(T)
    half = (p*np.sqrt(v)[:,None]/2).ravel()
    center = (q/2-m).ravel()
    N,A = N.ravel(),A.ravel()
    bounds = np.r_[0., edges, np.inf]
    out = np.zeros((2,len(bounds)-1))
    assert min(half)>0.
    for i in range(0,len(N),batch):
        s=slice(i,i+batch)
        x=np.clip((bounds[None,:]-center[s,None])/half[s,None],-1.,1.)
        mid=(x[:,1:]+x[:,:-1])*.5
        h=np.diff(x,axis=1)*.5
        avg=N[s,None]+A[s,None]*(3*mid**2-1+h**2)*.5
        n=2*h*avg
        k=2*h*(avg*np.maximum(0.,center[s,None]+half[s,None]*mid)
               +half[s,None]*A[s,None]*mid*h*h)
        out+=np.array([n.sum(axis=0),k.sum(axis=0)])
    return out


def interval_moments(lo, hi, width, rate, edges, batch=2048):
    """Integrate positive disjoint energy intervals including both tails."""
    bounds = np.r_[0., edges, np.inf]
    moments = np.zeros((2, len(bounds)-1))
    assert np.all(width > 0)
    for i in range(0, len(rate), batch):
        s = slice(i, i+batch)
        low = np.maximum(lo[s, None], bounds[:-1])
        high = np.minimum(hi[s, None], bounds[1:])
        n = ((rate[s]/width[s])[:, None])*np.maximum(0., high-low)
        # Inactive bins must not enter traces; inf is never part of an overlap.
        kmean = np.where(n > 0., (low+high)*.5, 0.)
        moments[0] += n.sum(axis=0)
        moments[1] += (n*kmean).sum(axis=0)
    return moments


class Visible:
    def __init__(self, order=(8,22,16), polar=24):
        self.direct = Source(order, dark_mass=m)
        self.heavy = Source(order, dark_mass=M)
        self.polar = polar

    def rates(self, T, q, a):
        q = np.asarray(q, dtype=float)
        assert len(q)>1 and min(q)>=0 and np.all(np.diff(q)>0)
        ke, kp = [kinetic(q/a, mass) for mass in (m, mp)]
        direct = direct_intervals(T, self.direct, ke)
        E, W = self.heavy.particles(T, n=self.polar)
        p = np.sqrt((E-M)*(E+M))
        width = 2*p*P/M
        pe_high = (Es*E+P*p)/M
        pe_low = (Es*Es+(mp*p/M)**2)/pe_high
        pe = interval_moments(pe_low-mp, pe_high-mp, width, W, kp)
        Ef = M-Es
        ce_high = (Ef*E+P*p)/M
        ce_low = (Ef*Ef+(m*p/M)**2)/ce_high
        cf = interval_moments(ce_low-m, ce_high-m, width, W, ke)
        count = W.sum()
        energy = W@E
        assert np.max(abs(np.array([pe[0].sum(), cf[0].sum()])/count-1)) < 3e-10
        assert abs((pe[1].sum()+cf[1].sum()+mp*pe[0].sum()+m*cf[0].sum())/energy-1) < 3e-10
        assert np.all(np.array([direct, cf, pe]) >= 0.)
        return np.array([direct, cf, pe])


def self_test(eos, T=.045):
    from empty_lepton_pair_spectrum_reference import reduced_moments
    from expanding_tagged_injection_transport_reference import prompt_source
    v = Visible()
    q = T*np.geomspace(.16, 16, 38)
    a = 1.04
    x = v.rates(T, q, a)
    mass = np.array([[m], [m], [mp]])
    n, e = x[:,0], x[:,1] + mass*x[:,0]
    nd, ed, _ = reduced_moments(T, v.direct)
    nh, eh, _ = reduced_moments(T, v.heavy)
    h = np.sqrt(eos.state(T)[0]/3)/MPL
    old = prompt_source(eos, T, np.r_[0.,q], a)*h*eos.state(T)[2]
    res = dict(direct_count=abs(n[0].sum()/(2*nd)-1),
               direct_energy=abs(e[0].sum()/ed-1),
               each_daughter_count=float(max(abs(n[1:].sum(axis=1)/(2*nh)-1))),
               joint_prompt_energy=abs(e[1:].sum()/eh-1),
               scalar_histogram_L1=float(np.sum(abs(old-n[2]))/(2*nh)))
    assert max(res.values()) < 3e-10
    return dict(T_GeV=T,a=a,q_over_T=[q[0]/T,q[-1]/T],defects=res,
                excluded_count_fraction=(n[:,[0,-1]]/n.sum(axis=1)[:,None]).tolist(),
                excluded_energy_fraction=(e[:,[0,-1]]/e.sum(axis=1)[:,None]).tolist())


def integrated_count_test(eos, raw, TR=.05, TF=.045):
    """No reactions; compare from empty to the independent source-only yield."""
    from scipy.optimize import brentq
    source = Visible(order=(9,24,16))
    q = TR*np.geomspace(.16,16,29)
    s0 = eos.state(TR)[2]
    extent = np.log(s0/eos.state(TF)[2])/3
    target = np.array([*raw["stable_yields"],raw["stable_yields"][1]])
    target_energy = target*np.array(raw["redshifted_means"])
    ms = np.array([m,m,mp])[:,None]
    final_a = np.exp(extent)
    estimates = []
    for n in (6, 12):
        dx = extent/n
        counts = np.zeros((3,3)); energy = np.zeros(3)
        for i in range(n):
            a = np.exp(dx*(i+.5))
            t = brentq(lambda t:eos.state(t)[2]/s0-a**-3,TF,TR)
            h = np.sqrt(eos.state(t)[0]/3)/MPL
            r = source.rates(t,q,a)
            nn = r[:,0]
            factor=dx*a**3/h
            counts += np.stack([nn[:,0], nn[:,1:-1].sum(axis=1), nn[:,-1]], axis=1)*factor
            # Nonlinear redshift need not preserve the one-node
            # approximation to the within-shell (and tail) mean energy.
            kin=np.divide(r[:,1],nn,out=np.zeros_like(nn),where=nn>0)
            p2=kin*(kin+2*ms)
            final_energy=np.sqrt(ms*ms+p2*(a/final_a)**2)
            energy+=factor*np.sum(nn*final_energy,axis=1)
        yield_split = counts/s0
        defect = max(abs(yield_split.sum(axis=1)/target-1))
        estimates.append(dict(steps=n,maximum_relative_yield_defect=float(defect),
                              approximate_redshifted_energy_yields=(energy/s0).tolist(),
                              approximate_maximum_relative_energy_yield_defect=float(
                                  max(abs((energy/s0)/target_energy-1))),
                              tails_and_kept_yields=yield_split.tolist()))
    assert estimates[1]["maximum_relative_yield_defect"] < min(4e-4,estimates[0]["maximum_relative_yield_defect"])
    assert estimates[1]["approximate_maximum_relative_energy_yield_defect"]<5e-4
    return estimates


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('eos')
    p.add_argument('--independent-spectra')
    a = p.parse_args()
    eos = StandardModel(a.eos)
    result = self_test(eos)
    if a.independent_spectra:
        with open(a.independent_spectra) as f:
            rows = json.load(f)['rows']
        raw = next(r for r in rows if r['TR']==.05 and abs(r['TF']-.045)<1e-8)
        result["source_only_from_empty_yield"] = integrated_count_test(eos,raw)
    print(json.dumps(result, indent=2, allow_nan=False))
