"""Transfer positive exclusive three-pion dispersion moments to the dark count ledger."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import M, m
from warm_babar_dispersive_gap_reference import MU
from warm_direct_vacuum_absolute_source_envelope_reference import (
    POINT, clock, log_slope_minimum, mpi, mu_kernel_lower, read, weight)
from warm_neutral_pion_photon_budget_reference import m0

SND_AND_BABAR=(
    ("SND",((.62,1.10,42.96,.06,.46),(1.10,1.975,2.99,.02,.08))),
    ("BABAR",((.62,1.10,42.91,.14,.55,.09),(1.10,2.,2.95,.03,.16))),
)
LOW=2*mpi+m0


def transfer(a,b,T,U,n,normalization_padding=1.):
    x=np.geomspace(a*a,b*b,n+1)
    l,h=x[:-1],x[1:]
    B=weight(l,h,T,U,False,mass=M)*2+weight(l,h,T,U,False,mass=m)
    den=mu_kernel_lower(h)/(4*np.pi**3)/normalization_padding
    v,w=roots_legendre(88)
    v=(v+1)/2
    check=((w/2)*(v*v*(1-v)))@(1/(v[None,:]**2+(1-v[None,:])*h[:,None]/MU**2)).T
    assert np.min(check/(den*4*np.pi**3))>=normalization_padding
    Fr=B*h/POINT/den
    k=np.argmax(Fr)
    return dict(energy_window_GeV=[a,b],gain_in_weighted_yield_per_dispersion_upper=float(Fr[k]),
                 controlling_s_interval_GeV2=[l[k],h[k]],lower_muon_kernel_relative_accuracy=float(np.min(den*4*np.pi**3/check)))


def flat(a,b,T,U,n):
    x=np.geomspace(a*a,b*b,n+1)
    W=lambda ma:weight(x[:-1],x[1:],T,U,False,mass=ma)
    return float(np.diff(x)@(2*W(M)+W(m)))


def run(eos,bracket,baseline):
    bg=Background(eos)
    lim=read(baseline)
    unit=read(bracket)["Y_to_Omega_h2"]
    check=[log_slope_minimum(bg.gs,3),log_slope_minimum(bg.gr,4)]
    out=[]
    for row in lim["rows"]:
        top=row["T_reheat_GeV"]
        fixed=row["controls"][-1]["total_count_scalar_equivalent_upper"]
        controls=[]
        for n in (240,480):
            T,U=clock(bg,3*n,top)
            U=U/((-np.expm1(-mpi/T))**2*(-np.expm1(-m0/T)))
            below=unit*flat(LOW,.62,T,U,n)
            variants=[]
            for label,measurements in SND_AND_BABAR:
                high_start=measurements[-1][1]
                above=unit*flat(high_start,4.,T,U,n)
                # SND writes a denominator 4*pi*alpha^2/s in (14)
                # with the conventional R prefactor. Allow a factor of
                # three without presuming it is a typesetting mistake.
                normalization_padding=3. if label=="SND" else 1.
                terms=[]
                for a,b,mean,*errors in measurements:
                    cap=(mean+2*sum(errors))*1e-10
                    f=transfer(a,b,T,U,n,normalization_padding)
                    q=f["gain_in_weighted_yield_per_dispersion_upper"]*cap*unit
                    terms.append(dict(**f,conditional_dispersion_moment_cap=cap,
                                      weighted_omega_upper=q))
                dq=sum(x["weighted_omega_upper"] for x in terms)
                total=fixed+dq
                assert .0<dq and total<.060
                variants.append(dict(experiment=label,bands=terms,normalization_padding=normalization_padding,
                        high_R1_interval_GeV=[high_start,4.],high_gap_weighted_omega_per_unit_R=above,weighted_omega_summed_measured_bands_upper=dq,
                        total_upper_including_only_measured_bands=total,
                        sub_062_uniform_R_budget_before_reaching_half=(.060-total)/below,
                        numerical_count_comparison_with_unresolved_low_and_high_R1_and_sqrt_s_cap_4GeV=total+below+above))
            controls.append(dict(invariant_panels=n,control_cases=variants,
                                 low_unmeasured_omega_ceiling_per_R=below))
        assert controls[1]["low_unmeasured_omega_ceiling_per_R"]<controls[0]["low_unmeasured_omega_ceiling_per_R"]
        for p,q in zip(controls[0]["control_cases"],controls[1]["control_cases"],strict=True):
            assert q["weighted_omega_summed_measured_bands_upper"]<p["weighted_omega_summed_measured_bands_upper"]
        out.append(dict(T_reheat_GeV=top,baseline_weighted_omega_upper=fixed,controls=controls))
    return dict(scope="Exclusive charged+charged+neutral three-pion contribution to the restricted weighted count, "
                      "not inclusive hadrons. Independent SND and BABAR alternatives use caps equal to the quoted "
                      "central values plus twice the error sums, not confidence intervals. SND carries an additional "
                      "factor-three safety multiplier for the normalization printed in its Eq. (14). The "
                      "optional low and high R=1 comparisons use disjoint windows, omit energy above "
                      "4 GeV, do not assert R=1 in those regions, and do not "
                      "include other species or finite-temperature spectral changes.",
                source_snd="https://arxiv.org/html/2603.01635#S7.T7",
                source_babar="https://doi.org/10.1103/PhysRevD.104.112003",
                threshold_GeV=LOW,clock_min_slopes=check,rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for f in ("eos","bracket","baseline"):
        p.add_argument(f,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.bracket,a.baseline),indent=2,allow_nan=False))
