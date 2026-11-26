"""Energy headroom and a restricted measured pion source on the warm clock."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import kve, roots_legendre

from kloe_pion_band_portal_reference import POINT_MUON, table
from late_radiation_background_reference import MPL, TabulatedRadiation
from three_species_decay_onsager_reference import M, m, mp
from measured_partial_hadron_reheat_fold_reference import integrands
from timelike_portal_onsager_coverage_reference import factor as spectral_prefactor
from warm_from_empty_recycling_envelope_reference import (
    S0, RHOC_H2, COUP2, VECTOR_MASS, VECTOR_GAMMA, energy_rate_coefficient,
    temperature_time_upper, upper_pair_Q,
)


def pion_positive_Q(T, data, splits, plus_error=True):
    """Nonnegative upper cells on the prescribed binwise vacuum source."""
    s,sigma,err = data
    b = (s-.005)[:,None] + .01*np.linspace(0.,1.,splits+1)
    l,h = b[:,:-1],b[:,1:]
    assert min(l.ravel())>.099 and max(h.ravel())<1.
    u = np.sqrt(l)/T
    k0,k1 = [kve(j,u)*np.exp(-u) for j in (0,1)]
    W = spectral_prefactor * (sigma+plus_error*err)[:,None] * h / POINT_MUON
    W *= h**1.5 / ((h-VECTOR_MASS**2)**2 + (VECTOR_MASS*VECTOR_GAMMA)**2)
    Q=N=0.
    for mass, prompt in ((m,2),(M,4)):
        p = W*(h+2*mass**2)*np.sqrt(1-4*mass**2/h)*(h-l)
        N += prompt*float(np.sum(p*(T*k1)))
        Q += float(np.sum(p*(2*T*T*k1+T*np.sqrt(h)*k0)))
    return np.array([Q,N])

def verify_pion_anchor(T,data,splits):
    s,cs,err=data
    x,w=roots_legendre(24)
    z=s[:,None]+.005*x
    W=(cs+err)[:,None]*z/POINT_MUON*(.005*w)
    input_row=(z.ravel(),W.ravel(),np.zeros(z.size))
    v=integrands(T,[input_row])[0]
    fixed=np.array([v[:,1].sum(),v[0,0]+2*v[1,0]])
    approx=pion_positive_Q(T,data,splits)
    q=approx/fixed
    assert np.all(q>=1.)
    return q.tolist()


def run_case(eos,K,measurement,n,du,splits,TR=.05,stop=500.,
             pion_rate_override=None,pion_tail_multiplier=None,heat_multiplier=1.,
             return_components=False):
    rad = TabulatedRadiation(eos,TR)
    lo = np.log(rad.state(stop)[0])
    hi = rad.initial
    ticks = rad.sg.x[(rad.sg.x>lo)&(rad.sg.x<hi)]
    u = np.unique(np.r_[np.linspace(lo,hi,n+1),ticks])[::-1]
    lgs,lgr = rad.sg(u),rad.rg(u)
    t = (rad.s0-3*u-lgs)/3
    Hph = np.pi/np.sqrt(45)*TR**2/MPL
    data = table(measurement)
    h = np.exp(.5*(np.minimum(lgr[:-1],lgr[1:])-np.log(2))+
               (2/3)*(lgs[0]-np.maximum(lgs[:-1],lgs[1:]))-1e-12)
    Ql = np.asarray([upper_pair_Q(T,du) for T in np.exp(u[:-1])])
    p_rate = pion_positive_Q if pion_rate_override is None else pion_rate_override
    Qpi = np.asarray([p_rate(T,data,splits) for T in np.exp(u[:-1])])
    assert np.all(np.diff(Qpi,axis=0)<=1e-200)
    assert np.all(np.diff(t)>0) and abs(t[-1]-np.log(stop))<1e-12
    # A convex analytic upper with fixed positive degree-of-freedom floors
    # at and after the cutoff; use a simple vacuum spectral maximum for
    # the far-cold pion source and the previous lepton massless ceiling.
    # R_eff <= max_j[(sigma_j + err_j) s_hi/POINT_MUON].
    Rmax=float(max((data[1]+data[2])*(data[0]+.005)/POINT_MUON))
    qcoeff = energy_rate_coefficient(.01045)[0]
    pion_qcoeff = Rmax*qcoeff/3
    if pion_tail_multiplier is not None:
        pion_qcoeff *= pion_tail_multiplier
    inv=(np.exp((lgs[0]-np.log(2))/3+1e-12)*TR)**9/(4*Hph*stop**4)
    l_tail,p_tail = qcoeff*inv,pion_qcoeff*inv
    dt = np.diff(np.exp(5*t))/(5*Hph*h)
    E = np.cumsum(dt*Ql[:,0])
    Epi = np.cumsum(dt*Qpi[:,0])
    El = E[-1]+l_tail
    Ep = Epi[-1]+p_tail
    Nl = dt@Ql[:,1] + l_tail/m
    Np = dt@Qpi[:,1] + p_tail/m
    W = np.diff(-np.exp(-t))/h
    w0 = W.sum()+1/stop
    w1 = W@E+El/stop
    w2 = W@(E**2)+El**2/stop
    heat = 288*COUP2*temperature_time_upper(eos,TR,stop)[1]/(np.pi**3*VECTOR_MASS**4)
    c = K*np.exp(2*heat*heat_multiplier)/(M*M*Hph)
    s_init = 2*np.pi**2/45*np.exp(lgs[0])*TR**3
    nrm = mp*S0/(RHOC_H2*s_init)
    # Permit unknown positive energy through this explicit budget only.
    # An early deposit is conservative for recycling on every later slice.
    def upper(e,with_count):
        return nrm*(Nl+with_count*e/m+c*(w2+2*e*w1+e*e*w0))
    headroom = brentq(lambda r: upper(r*El,True)-.060,0.,100.,
                      xtol=1e-13)*El
    pure_heat_headroom=brentq(lambda r: upper(r*El,False)-.060,0.,100.,
                             xtol=1e-13)*El
    et = El+Ep
    xx = E+Epi
    w1_phased = W@xx+et/stop
    w2_phased = W@(xx**2)+et**2/stop
    def with_pion_extra(e):
        return nrm*(Nl+Np+e/m+c*(w2_phased+2*e*w1_phased+e*e*w0))
    phased=with_pion_extra(0.)
    leftover=brentq(lambda r: with_pion_extra(r*El)-.060,0.,100.,
                     xtol=1e-13)*El
    inflated=upper(Ep,True)
    assert phased <= inflated and upper(0.,True) < phased
    assert abs(upper(headroom,True)-.060)<2e-14
    # Unlike observational statistics, the quoted per-bin +1 error is
    # a model input; this comparison makes no joint confidence claim.
    result=dict(panels=len(t)-1,du=du,pion_subcells=splits,TR_GeV=TR,
                model_lepton_energy_upper_GeV4=float(El),
                prescribed_pion_energy_upper_GeV4=float(Ep),
                prescribed_pion_energy_over_lepton=float(Ep/El),
                prescribed_pion_count_over_lepton=float(Np/Nl),
                pion_cell_upper_over_independent_anchor_fold=verify_pion_anchor(TR,data,splits),
                lepton_only_scalar_Omega_ceiling=float(upper(0.,True)),
                phased_measured_pion_and_lepton_scalar_Omega_ceiling=float(phased),
                all_pion_energy_pushed_earliest_scalar_Omega_ceiling=float(inflated),
                extra_portal_energy_for_ceiling_to_touch_half_GeV4=float(headroom),
                extra_portal_energy_over_lepton_to_touch_half=float(headroom/El),
                extra_energy_after_including_pions_to_touch_half_GeV4=float(leftover),
                extra_energy_after_including_pions_over_lepton=float(leftover/El),
                pure_heat_energy_over_lepton_to_touch_half=float(pure_heat_headroom/El))
    if return_components:
        result["baseline_components"]=dict(direct_Omega=float(nrm*(Nl+Np)),
                 recycle_Omega=float(nrm*c*w2_phased), heat_exponent=float(heat*heat_multiplier))
    return result


def run(eos,measurement,kernel):
    a = json.loads(Path(kernel).read_bytes())
    K = a["simultaneous_energy_density_coefficient_GeVminus2"]
    A=[run_case(eos,K,measurement,*j) for j in ((400,.02,4),(850,.009,8),(1200,.006,12))]
    assert A[-1]["lepton_only_scalar_Omega_ceiling"]<A[0]["lepton_only_scalar_Omega_ceiling"]
    return dict(scope="Initially empty 50 MeV instantaneous vacuum observable and "
                      "central EOS. KLOE (central plus per-bin quoted error) on the "
                      "prescribed piecewise-constant exclusive 0.10–0.95 GeV^2 band "
                      "only. Not an envelope of unknown hadronic rates, systematics, "
                      "thermal spectral changes or scalar external injection.",
                input_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest()
                              for p in (eos,measurement,kernel)],cases=A)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","measurement","kernel"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.measurement,a.kernel),indent=2,allow_nan=False))
