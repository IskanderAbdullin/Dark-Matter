"""Distance from the saturated thermal normal branch at fixed fast light count and dark energy."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit,roots_legendre,zeta
from scipy.optimize import brentq

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import m,mp,M
from warm_prompt_reheating_expanding_slow_reference import visible


def measure(T,n):
    u,w=roots_legendre(n)
    # Include an explicit infrared panel and resolve m_phi/T. The
    # saturated mu=m scalar occupation is integrated without a floor.
    J=np.unique(np.r_[0.,np.geomspace(1e-8,82.,36),mp/T])
    q=T*((J[:-1,None]+J[1:,None])/2+np.diff(J)[:,None]*u/2)
    mass=np.array((m,mp,M))[:,None,None]
    E=np.hypot(q[None,:,:],mass)
    weight=np.array([4,1,4])[:,None,None]*w*np.diff(J)[:,None]*T*q*q/(4*np.pi**2)
    phi=1/np.expm1(q*q/((E[1]+mp)*T))
    assert min(phi.ravel())>=0
    return E,weight,phi


def state(T,L,n):
    E,w,F=measure(T,n)
    shift=mp/T
    def batches(eta):
        return np.array([expit(eta-E[0]/T),F,expit(eta+shift-E[2]/T)])
    def diag(eta):
        z=w*batches(eta)
        return z.sum(axis=(1,2)),np.sum(z*E)
    eta=brentq(lambda x:diag(x)[0][[0,2]].sum()/L-1.,-70.,max(20.,m/T+50.),xtol=1e-13)
    a,b=diag(eta)
    return a,b,eta


def compare(T,L,S,U,orders=(22,38)):
    out=[]
    for n in orders:
        f=lambda x:state(x,L,n)[1]/U-1.
        low=T/2
        while f(low)>0.:
            low/=2
            assert low>T*1e-4
        Tc=brentq(f,low,T,xtol=1e-13)
        q,E,eta=state(Tc,L,n)
        assert abs(E/U-1)<2e-11 and abs((q[0]+q[2])/L-1)<2e-11 and Tc<T
        cap=q[1]+q[2]
        out.append(dict(order=n,T_saturated=Tc,scalar_invariant_at_saturation=cap,
                        scalar_invariant_over_capacity=S/cap,cap_heavy_fraction=q[2]/q.sum(),eta_light_saturated=eta))
    assert abs(out[0]["scalar_invariant_over_capacity"]/out[-1]["scalar_invariant_over_capacity"]-1)<1e-8
    return out


def independent_massless_limit():
    T=.017
    E,w,f=measure(T,40)
    N=(w[1]*f).sum()
    rho=(w[1]*f*E[1]).sum()
    # At mu=m the leading small m/T changes are the chemical-potential
    # derivatives of the massless integrals. Do not use this to fit
    # the finite-mass phase boundary.
    predicted=(mp/T)*np.array([zeta(2)/zeta(3),zeta(3)/zeta(4)])
    observed=np.array([N/(zeta(3)*T**3/np.pi**2)-1,rho/(np.pi**2/30*T**4)-1])
    assert max(abs(observed/predicted-1))<1e-3
    return dict(T=T,fractional_number_and_energy_shifts=observed.tolist(),
                independent_small_mass_leading_terms=predicted.tolist())


def calc(eos,file):
    bg=Background(eos)
    x=json.loads(file.read_text())
    if "arms" in x:
        root=x["arms"][0]["rows"]
        data=[dict(T=p["T_dark"],visible_T=p["T_visible"],pop=p["physical_yields"],
                   inv=p["invariant_yields"],peak=p["peak"]) for p in root]
    else:
        root=x["history"][-1]["samples"]
        data=[dict(T=p["dark_T"],visible_T=p["visible_T"],pop=p["actual_counts"],inv=p["invariants"],peak=p["phase_peak"]) for p in root]
    rows=[]
    for i,d in enumerate(data):
        if i not in (0,len(data)//2,len(data)-1):
            continue
        T,v,p,inv,peak=(d[k] for k in ("T","visible_T","pop","inv","peak"))
        _,s=visible(bg,v)
        rho,weights,_=measure(T,38)
        m0=np.array([m,mp,M])
        eta=np.log(peak)+m0/T
        F=np.array([expit(eta[0]-rho[0]/T),1/np.expm1(rho[1]/T-eta[1]),expit(eta[2]-rho[2]/T)])
        pop=weights*F
        assert np.max(abs(pop.sum(axis=(1,2))/s/p-1))<2e-8
        U=np.sum(pop*rho)
        result=compare(T,inv[0]*s,inv[1]*s,U)
        rows.append(dict(T_visible=v,T_dark=T,peak=peak,controls=result))
    return dict(independent_massless_limit=independent_massless_limit(),scope="Compare to a saturated thermal normal phase at fixed light-equivalent number and total dark energy, "
                      "with the same fast heavy decay affinity. Does not determine kinetic condensation, the validity of "
                      "the thermal closure, or whether an evolving history reaches saturation.",rows=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("path",type=Path)
    a=p.parse_args()
    print(json.dumps(calc(a.eos,a.path),indent=2,allow_nan=False))
