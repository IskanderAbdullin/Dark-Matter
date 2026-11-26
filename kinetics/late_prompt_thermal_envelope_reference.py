"""Full-support adiabatic parent envelope for specified thermal FD/BE baths.

This is a bound for the supplied thermal distributions, not for arbitrary
high-energy nonthermal tails at the same number and energy moments.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import kve,logsumexp

from continuum_dark_normal_capacity_reference import fugacity,moments,rule
from three_species_decay_onsager_reference import m,mp,M,Es,P,width


def envelope(Tp,Tc,ep,ec):
    assert min(Tp,Tc)>0. and ep<0.
    th=max(Tp,Tc)
    z=M/th
    const=np.log(4*width*M*M*th/(2*np.pi**2))+ep+ec-np.log(-np.expm1(ep))
    th_diff=-(M-m-mp)/th
    return np.array([const+th_diff+np.log(kve(1,z)),
                     const+th_diff+np.log(M*kve(2,z))])


def quadrature(Tp,Tc,ep,ec,order=24):
    """Separate angular inverse-fold estimate and positive heat envelope."""
    T=max(Tp,Tc)
    x,w=rule(order)
    limits=np.array([0.,.3,1.5,5.,14.,34.,65.])
    e=M+T*((limits[:-1,None]+limits[1:,None])/2+np.diff(limits)[:,None]*x/2).ravel()
    weight=(T*np.diff(limits)[:,None]*w/2).ravel()
    p=np.sqrt((e-M)*(e+M))
    high=(e*Es+p*P)/M
    low=(Es*Es+(mp*p/M)**2)/high
    angular=low[:,None]+(high-low)[:,None]*(x+1)/2
    phi=(angular-mp)/Tp-ep
    chi=(e[:,None]-angular-m)/Tc-ec
    logf=-phi-np.log(-np.expm1(-phi))-np.logaddexp(0.,chi)
    logA=logsumexp(logf+np.log(w/2),axis=1)
    pref=np.log(4*M*width/(2*np.pi**2)*p*weight)
    return np.array([logsumexp(pref+logA),
                     logsumexp(pref+logA+np.log(high-low))])


def audit(initial_counts,rows):
    n0=np.array(initial_counts)
    out=[]
    for r in rows:
        if r["b"]<1.5 or r["b"]>26.: continue
        b=r["b"]
        n=n0/b**3
        Tp,Tc=r["equivalent_T_GeV"]
        ep=fugacity(Tp,mp,n[0],True,20)
        ec=fugacity(Tc,m,n[1],False,20)
        H=r["H_GeV"]
        K=moments(Tp,mp,ep,True)[1]
        common_bound=envelope(Tp,Tc,ep,ec)
        direct=quadrature(Tp,Tc,ep,ec)
        fine=quadrature(Tp,Tc,ep,ec,40)
        assert np.all(direct<=common_bound+1e-10) and np.all(fine<=common_bound+1e-10)
        norm=np.log([n[0]*H,K*H])
        out.append(dict(b=b,T_phi_GeV=Tp,T_fermion_GeV=Tc,
            log10_inverse_events_per_scalar_per_Hubble_envelope=float((common_bound[0]-norm[0])/np.log(10.)),
            log10_abs_heat_over_scalar_kinetic_Hubble_envelope=float((common_bound[1]-norm[1])/np.log(10.)),
            log10_angular_inverse_rate_per_scalar_per_Hubble=float((fine[0]-norm[0])/np.log(10.)),
            log10_angular_abs_heat_envelope_over_scalar_kinetic_Hubble=float((fine[1]-norm[1])/np.log(10.)),
            log10_quadrature_change=float(max(abs(fine-direct))/np.log(10.))))
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("capacity");p.add_argument("late_trace")
    a=p.parse_args()
    n=json.loads(Path(a.capacity).read_text())['rows'][-1]['species_counts_GeV3']
    r=json.loads(Path(a.late_trace).read_text())["steps"]
    out=audit(n,r)
    print(json.dumps(dict(rows=out,domain="specified full continuum thermal distributions; not a constraint on unknown nonthermal tails"),indent=2,allow_nan=False))
