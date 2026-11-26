"""Upper envelope retaining the specified late charged-lepton portal sources."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gamma, gammaincc, kve
from scipy.integrate import quad

from late_radiation_background_reference import TabulatedRadiation, MPL
from three_species_decay_onsager_reference import M, m, mp


VECTOR_MASS = 10.
VECTOR_GAMMA = .005305164626087453
COUP2 = .1**2 * 1e-5**2 * 4*np.pi/137.035999084
LEPTONS = (.00051099895, .1056583755, 1.77686)


def energy_rate_coefficient(T_cap, split_sqrt=1.):
    """Q_MB <= coefficient * T**9 for all T<=T_cap, in GeV."""
    # The fermion mass/phase factors obey (1+2r)*sqrt(1-4r)<=1.
    u = split_sqrt / T_cap
    assert u > 10. and split_sqrt < VECTOR_MASS
    integral = 2*2**5*gamma(4)*gamma(3) + 2**6*gamma(4)**2
    below = integral/(VECTOR_MASS**2-split_sqrt**2)**2
    # For nu=0,1, cosh(nu*x) <= exp(x*x/2), cosh(x)>=1+x*x/2.
    # Thus K_nu(u)<=e^-u sqrt(pi/(2(u-1))). Use an integer tail sum.
    tail = math.sqrt(np.pi/(2*(u-1.))) * (
        2*gamma(7)*gammaincc(7,u) + gamma(8)*gammaincc(8,u))
    BW_floor = (VECTOR_MASS*VECTOR_GAMMA)**2
    above = tail/BW_floor
    return 6*COUP2/(48*np.pi**5)*(below+above), above/below


def temperature_time_upper(eos,T,b_stop=500.):
    R = TabulatedRadiation(eos,T)
    Hph = np.pi/np.sqrt(45)*T*T/MPL
    lo = np.log(R.state(b_stop)[0])
    x = R.sg.x
    u = np.r_[lo,x[(x>lo)&(x<R.initial)],R.initial]
    sg,rg = R.sg(u),R.rg(u)
    sg_ref = sg[-1]
    gslo = np.minimum(sg[:-1],sg[1:])
    gshi = np.maximum(sg[:-1],sg[1:])
    rglo = np.minimum(rg[:-1],rg[1:])
    t = (R.s0-3*u-sg)/3
    assert np.all(np.diff(t)<0.) and abs(t[0]-np.log(b_stop))<1e-12
    # Pad the endpoint-based bounds outward, as in the isolated clock.
    T_scale = np.exp((sg_ref-gslo)/3+1e-12)
    H_scale = np.exp(.5*(rglo-np.log(2))+(2./3)*(sg_ref-gshi)-1e-12)
    assert min(H_scale)>1.
    after = np.exp((sg_ref-np.log(2.))/3+1e-12)
    # At b>=b_stop condition on H>=Hph/b^2 and at least two photon
    # entropy degrees of freedom on the supplied constant-entropy history.
    pow_terms = []
    for power,power_b in ((9,4),(5,3)):
        early = float(np.sum(np.diff(np.exp(-power_b*t))/power_b
                            * T_scale**power/H_scale))
        late = after**power / (power_b*b_stop**power_b)
        pow_terms.append(T**power/Hph*(early+late))
    return pow_terms


def direct_MB(T):
    """Independent nonzero-mass invariant fold; diagnostic, not an enclosure."""
    def pair(mass,lepton):
        u0 = 2*max(mass,lepton)/T
        f = lambda z: (1+2*z)*np.sqrt(max(0.,1-4*z))
        def kernel(x):
            u=u0+x
            s=(T*u)**2
            return (u**6*f(mass**2/s)*f(lepton**2/s)*np.exp(-x)
                    *(2*kve(1,u)+u*kve(0,u))/
                    ((s-VECTOR_MASS**2)**2+(VECTOR_MASS*VECTOR_GAMMA)**2))
        return np.exp(-u0)*quad(kernel,0.,80.,epsabs=1e-28,epsrel=2e-9)[0]
    return COUP2*T**9/(48*np.pi**5)*sum(pair(d,l) for d in (m,M) for l in LEPTONS)


def run(eos,isolated_path,exterior_path):
    A=json.loads(Path(isolated_path).read_bytes())
    B=json.loads(Path(exterior_path).read_bytes())
    kappa=288*COUP2/(np.pi**3*VECTOR_MASS**4)
    ans=[]
    for s,t in zip(A["cases"],B["cases"],strict=True):
        T,rho,N = s["T_GeV"],s["represented_energy_GeV4"],s["represented_number_GeV3"]
        assert s["q"]==t["q"] and abs(T/t["anchor_T_GeV"]-1.)<1e-12
        Y=s["conservative_initial_Y"]
        goal=A["target_total"]
        rho_complete=rho*(1+t["bounds"][-1]["added_energy_over_core"])
        int_pair,int_elastic = temperature_time_upper(eos,T)
        QC,tail_over_contact=energy_rate_coefficient(T)
        j_pair=QC*int_pair
        j_heat=kappa*int_elastic
        # Positive pressure and no other positive external energy: use
        # the Gronwall bound for all times; grants all future pair input early.
        rho_cap=np.exp(j_heat)*(rho_complete+j_pair)
        # Each visible pair supplies at most one final count per m of energy,
        # including both daughters per directly injected heavy particle.
        direct_n=j_pair/m
        assert M>=2*m
        time=s["bounds"][1]["clock"]
        recycle_n=(A["combined_K_GeVminus2"]/M**2*rho_cap**2
                   *time["upper_integral_times_Hph"]/time["photon_anchor_GeV"])
        mul=1+(direct_n+recycle_n)/N
        X=(Y*mul-mp/m*goal)/(1-mp/m)
        cmp=direct_MB(T)/(QC*T**9)
        assert 0<cmp<1.
        ans.append(dict(q=s["q"],anchor_T_GeV=T,pair_energy_upper_over_core=j_pair/rho,
                        elastic_integrating_factor_exponent=j_heat,
                        direct_count_upper_over_core=direct_n/N,
                        recycling_count_upper_over_core=recycle_n/N,
                        all_particles_multiplier_upper=mul,
                        scalar_Omega_h2_upper_at_reference_total=X,
                        scalar_share_upper_at_reference_total=X/goal,
                        massless_pair_bound_resonant_tail_fraction=tail_over_contact,
                        direct_massive_MB_anchor_under_upper=cmp,
                        additional_early_rho_budget_GeV4=rho_cap-rho_complete))
    return dict(scope="Known cohorts and arbitrary redistribution, positive late vacuum "
                      "charged-lepton pair injection and elastic heating, specified isolated "
                      "internal network, fixed visible entropy and expansion floor. "
                      "No late hadronic sources, new dark channels, or untracked external inflow.",
                target=A["target_total"],kappa_over_T5=kappa,cases=ans)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","isolated","exterior"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.isolated,a.exterior),indent=2,allow_nan=False))
