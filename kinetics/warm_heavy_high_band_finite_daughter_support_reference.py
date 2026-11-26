"""Conditional high heavy tail bound when both light daughter supports stay finite."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.special import gamma,gammaincc

from late_radiation_background_reference import TabulatedRadiation
from three_species_decay_onsager_reference import M,m,mp
from warm_archived_bose_pair_weighting_reference import snapshot
from warm_fermion_survival_with_leptons_reference import subpole
from warm_late_visible_recycling_envelope_reference import COUP2,VECTOR_MASS,VECTOR_GAMMA


def read(p):
    x=Path(p).read_bytes()
    return json.loads(gzip.decompress(x) if str(p).endswith(".gz") else x)


def max_produced(P,a,b):
    u,v=np.hypot(P,a),np.hypot(P,b)
    # Treat all collisions as if the tagged heavy is the only output.
    # As a function of b_exp^2, this upper squared comoving momentum
    # is concave. Its derivative is already negative at the anchor.
    derivative=a*a+b*b-M*M+(a*a*v*v+b*b*u*u)/(u*v)
    assert derivative<0
    return float(np.sqrt(max(0.,(u+v)**2-M*M)))


def high_visible(T,temp_factor,I3,I4,cut,Kpi):
    # Integral of E exp(-E/T) over p >= cut/b, using E <= M+p.
    def with_scale(a):
        t=a*temp_factor*T
        z=cut/t
        tail3=gamma(3)*gammaincc(3,z)
        tail4=gamma(4)*gammaincc(4,z)
        return (4/np.pi**4)*(M*t**6*tail3*I3+t**7*tail4*I4)
    V=VECTOR_MASS
    w=V*VECTOR_GAMMA
    s=np.hypot(V*V,w)
    vmax=3*COUP2*s/(4*np.pi*((s-V*V)**2+w*w))
    low=(subpole(1.)+Kpi)*with_scale(1.)
    # Above sqrt(s)=1, reserve half the thermal exponent for that
    # invariant and integrate the remaining Maxwell weights at 2T.
    high=vmax*np.exp(-1/(2*T))*with_scale(2.)
    return low,high


def run(eos,pions,isolated,aggregate,light,shapes):
    A,I,J,L=map(read,(pions,isolated,aggregate,light))
    out=[]
    for a,i,j,l,p in zip(A["cases"],I["cases"],J["cases"],L["cases"],shapes,strict=True):
        assert a["q"]==i["q"]==j["q"]==l["q"]
        S=snapshot(p)
        P=S["lab_momentum_edges_GeV"][-1]; T=S["T_visible_GeV"]
        assert abs(T/i["T_GeV"]-1)<1e-12
        row=a["controls"][-1]
        E=row["new_available_energy_upper_GeV4"]; Kpi=row["pion_pair_sigma_v_cap_with_stated_support"]
        R=TabulatedRadiation(eos,T)
        aT=(np.exp(float(R.sg(R.initial)))/2)**(1/3)
        clock=i["bounds"][1]["clock"]; stop=clock["b_stop"]
        hp,h=clock["photon_anchor_GeV"],clock["floor_min"]
        I3=clock["upper_integral_times_Hph"]/hp
        I4=((1-stop**-2)/h+stop**-2)/(2*hp)
        maxima=[max_produced(P,*z) for z in ((m,mp),(m,m),(mp,mp))]
        pj=next(z for z in j["rows"] if z["light_cap"]==.5)
        sub=next(z for z in pj["controls"] if z["cut_GeV"]==1.)
        assert max(maxima)<1.
        v,w=high_visible(T,aT,I3,I4,1.,Kpi)
        frac=(v+w)/E
        total=sub["strongest_bound_for_below_cut"]+frac
        assert total<l["diagnostic_heavy_energy_fraction_for_direct_dose_half"]
        out.append(dict(q=a["q"],assumed_equal_comoving_daughter_cap_float_GeV=P,
                         internal_produced_heavy_cap_GeV=maxima,
                         high_cut_GeV=1.,
                         high_band_visible_energy_GeV4=[v,w],high_fraction_upper=frac,
                         combined_heavy_energy_fraction_upper=total,
                         light_direct_diagnostic_limit=l["diagnostic_heavy_energy_fraction_for_direct_dose_half"]))
    return dict(scope="Both dark daughter comoving supports assumed capped by the "
                      "archive edges throughout this conditional comparison, all heavy "
                      "modes initially empty, no heavy elastic exchange. Actual new or "
                      "untracked exterior daughter populations are not included. High "
                      "visible heavy injection is bounded, not set to zero. This is not "
                      "a proof of the assumed support or actual light occupation cap.",
                cases=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("eos","pions","isolated","aggregate","light"):
        p.add_argument(x)
    p.add_argument("shapes",nargs=2)
    a=p.parse_args()
    print(json.dumps(run(*vars(a).values()),indent=2,allow_nan=False))
