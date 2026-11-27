"""Conditional late stable floor without fixed momenta or a collisionless assumption."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from archived_scalar_pair_reaction_reference import amplitude
from conditional_relic_normalization_reference import conversion
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp,M,y
from warm_prompt_reheating_expanding_slow_reference import visible
from warm_late_visible_recycling_envelope_reference import VECTOR_MASS, VECTOR_GAMMA, COUP2


def run(eos,prior):
    bg=Background(eos)
    d=json.loads(prior.read_text())
    x,w=roots_legendre(88)
    s=np.geomspace(4*m*m*(1+1e-6),100.,500)
    # One opposite-charge event, identical final scalar factor included.
    direct=(amplitude(s[:,None],x)@w)*np.sqrt(1-4*mp*mp/s)/(32*np.pi*s)
    cap=y**4/(8*np.pi*M*M)
    assert max(direct)<cap
    p=VECTOR_MASS**2
    w2=(VECTOR_MASS*VECTOR_GAMMA)**2
    star=np.sqrt(p*p+w2)
    maxgeom=(star+p)/(2*w2)
    vector=3*COUP2*maxgeom/(4*np.pi)
    assert abs(star/((star-p)**2+w2)/maxgeom-1)<2e-15
    om=conversion(eos)["Omega_h2_per_GeV_times_yield"]*m
    refs=[]
    for st in (d["arms"][-1]["trajectory"][0],d["arms"][-1]["trajectory"][-1]):
        tv=st["T_visible"]
        Y=st["yields"][0]
        _,sv=visible(bg,tv)
        Hph=np.sqrt(np.pi**2/45)*tv*tv/MPL
        C=cap*sv*Y/(2*Hph)
        V=vector*sv*Y/(2*Hph)
        raw=Y*om
        # The late inverse rate, if present, is positive in the light
        # count equation and does not weaken this differential floor.
        F=np.sqrt((raw/.120-1)/C)-1
        Fv=np.sqrt((raw/.120-1-V)/C)-1
        entries=[]
        for f in (.5,1.,4.,10.,30.,100.):
            k=1/(1+C*(1+f)**2)
            entries.append(dict(scalar_occupation_cap=f,retained_light_lower_fraction=k,
                                source_referenced_Omega_h2_lower=raw*k,
                                with_three_leptons_width_uniform_Omega_h2_lower=raw/(1+V+C*(1+f)**2)))
        refs.append(dict(T_visible_GeV=tv,initial_Y_light=Y,initial_Omega_h2_if_counts_freeze=raw,
                         exposure_coefficient_per_Bose_multiplier=C,uniform_vector_exposure_coefficient=V,
                         necessary_Bose_cap_for_floor_to_drop_to_point12=float(F),
                         with_three_leptons_necessary_Bose_cap_for_bound_to_drop_to_point12=float(Fv),controls=entries))
    return dict(scope="Only number-changing route for stable light fermions is light-antilight to scalar pair "
                      "and its inverse, optionally also to and from three charged leptons with a fixed vacuum vector width. Elastic scattering and particle momenta need not be bounded, "
                      "but the scalar occupation has an assumed uniform upper cap at all later times. "
                      "Photon-only expansion floor assumes T_gamma>=anchor/b; the chosen normalizer "
                      "uses conserved visible entropy. A prior or later nonstandard cosmology, extra loss "
                      "channels, or unbounded occupations are not excluded or solved here.",
                invariant_vacuum_upper_GeVminus2=cap,vector_width_uniform_upper_GeVminus2=vector,sampled_max_over_analytic_upper=float(max(direct)/cap),anchors=refs)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("eos",type=Path);p.add_argument("path",type=Path)
    a=p.parse_args();print(json.dumps(run(a.eos,a.path),indent=2,allow_nan=False))
