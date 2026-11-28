"""A pointwise allowance for accumulating light-fermion pair source."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp,y
from warm_fermion_survival_with_leptons_reference import subpole
from warm_late_visible_recycling_envelope_reference import COUP2,VECTOR_MASS,VECTOR_GAMMA


def read(p):
    return json.loads(Path(p).read_bytes())


def direct_dose(eos,T,E,Kpi,weight,margin=.5):
    r=TabulatedRadiation(eos,T)
    a=(np.exp(float(r.sg(r.initial)))/2)**(1/3)
    # A tagged light final spin sums the two partner spins. Dropping
    # the initial identical-scalar factor makes label=2 conservative.
    label=2
    beta1=np.sqrt(1-(2*mp/m)**2)
    beta2=np.sqrt(1-(mp/(2*m))**2)
    hardsoft=(2*np.pi**2*E)**2/m**4*(16/beta1**2+4/beta2)
    Sscalar=(label/(16*np.pi**3))*2*y**4/M**2*hardsoft/m * weight(6)
    # For the reverse heavy -> light conversion, permit f_heavy <=
    # half per spin up to the first exit from the light cap (the
    # independent heavy barrier must hold up to that event).
    def from_heavy(r):
        n=r*E/M
        J=.5*(3*np.pi**2*n)**(2/3)*.5**(1/3)
        d=4*((M+m)**2-mp*mp)/(M*M-m*m)**2
        hi=(J/M)**2+(d*(M*J+np.pi**2*n))**2
        return (32/(16*np.pi**3))*y**4/2*hi/m*weight(4)
    Sheavy=from_heavy(1.)
    # Bound visible direct production at a tagged light mode by
    # reversing the kernel with a Boltzmann partner density.
    # Below the split, enforce true support of the pion source.
    def Svis(cut):
        assert cut>=1.
        nb=2*(a*T)**3/np.pi**2
        v=VECTOR_MASS; w=v*VECTOR_GAMMA
        s=np.hypot(v*v,w)
        vmax=3*COUP2*s/(4*np.pi*((s-v*v)**2+w*w))
        lo=nb*np.exp(-m/T)*(subpole(cut)+Kpi)*weight(3)
        # Keep half the Boltzmann factor to integrate the other
        # outgoing fermion. The other half suppresses the invariant.
        high=8*nb*np.exp(-m/(2*T)-cut/(2*T))*vmax*weight(3)
        return lo+high
    cuts=np.linspace(1.,1.75,60)
    C=np.array([Svis(t) for t in cuts])
    v=float(min(C))
    total=float(Sscalar+Sheavy+v)
    assert 0<Sscalar+v<margin<total
    critical=brentq(lambda f:Sscalar+from_heavy(f)+v-margin,0.,1.,xtol=1e-14)
    controls=[dict(heavy_energy_fraction_cap=r,pair_dose=float(Sscalar+from_heavy(r)+v))
              for r in (0.,.01,.02,.04,.05,.10,1.)]
    return dict(total_direct_gain_dose_upper=total,
                      scalar_pair_contribution=float(Sscalar),heavy_pair_contribution=float(Sheavy),
                      visible_pairs_contribution=v,chosen_visible_split_GeV=float(cuts[np.argmin(C)]),
                      diagnostic_heavy_energy_fraction_for_direct_dose_half=critical,
                      test_energy_fraction_caps=controls,
                      residual_margin_to_half_ignoring_initial_and_other_sources=margin-total)


def run(eos,snapshot,pions,isolated):
    A,P,I=map(read,(snapshot,pions,isolated))
    cases=[]
    for f,p,i in zip(A["cases"],P["cases"],I["cases"],strict=True):
        q=f["q"]
        assert p["q"]==i["q"]==q
        E=p["controls"][-1]["new_available_energy_upper_GeV4"]
        Kpi=p["controls"][-1]["pion_pair_sigma_v_cap_with_stated_support"]
        clock=i["bounds"][1]["clock"]
        T=f["anchor_T_GeV"]; h=clock["floor_min"]; hp=clock["photon_anchor_GeV"]; stop=clock["b_stop"]
        def weight(power):
            return ((1-stop**(2-power))/h+stop**(2-power))/(hp*(power-2))
        cases.append(dict(q=q,**direct_dose(eos,T,E,Kpi,weight)))
    return dict(scope="Partial pointwise light-occupation upper gain after 10.44 MeV. "
                      "Inputs assume the per-spin heavy cap and use "
                      "the stated total-energy/visible clock. Direct pair sources only. "
                      "Does not include gains from the one-to-two decay/inverse network "
                      "or elastic redistribution. Not a proof of the light cap.",
                cases=cases)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("eos","snapshot","pions","isolated"):
        p.add_argument(k)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.snapshot,a.pions,a.isolated),indent=2,allow_nan=False))
