"""Exact stable coproduction ledger and restricted two-to-two count headroom."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from archived_collisionless_streaming_reference import Background
from empty_lepton_pair_spectrum_reference import Source
from sm_dark_entropy_exchange_reference import MPL
from three_species_decay_onsager_reference import m,mp
from warm_light_occupation_stage_audit_reference import read


def emit_light(eos,tops):
    b=Background(eos)
    rows=[]
    for rule,nt in (((9,28,24),22),((14,46,32),36)):
        src=Source(rule,tail=75.,dark_mass=m)
        w,x=roots_legendre(nt)
        yield_=[]
        for l,r in zip(np.log([.001,*tops[:-1]]),np.log(tops),strict=True):
            acc=0.
            # Resolve the kinks in EOS derivatives before integrating.
            splits=np.r_[l,b.u[(b.u>l)&(b.u<r)],r]
            for a,c in zip(splits[:-1],splits[1:],strict=True):
                for u,du in zip((c+a)/2+(c-a)*w/2,(c-a)*x/2,strict=True):
                    T=np.exp(u)
                    S=(2*np.pi**2/45)*T**5*np.exp(b.gs(u)+.5*b.gr(u))*np.sqrt(np.pi**2/90)/MPL
                    acc+=du*(1+b.gs.derivative()(u)/3)*2*src.event_moments(T)[3].sum()/S
            yield_.append(acc)
        rows.append(np.cumsum(yield_))
    a,z=rows
    assert np.max(abs(a/z-1.))<2e-5
    return z


def run(eos,bracket):
    r=read(bracket)
    data=r["rows"]
    tops=[x["T_reheat_GeV"] for x in data]
    light=emit_light(eos,tops)*r["Y_to_Omega_h2"]
    ratio=m/mp
    total=.12
    out=[]
    for x,y in zip(data,light,strict=True):
        s=sum(x["central_lepton_and_KLOE_Omega_h2"])
        stable_parent=s*ratio
        stable_light=y*ratio
        stable=stable_parent+stable_light
        fraction=(total-s)/stable
        assert 0<fraction<1
        upper=x["controls"][-1]["upper"]
        out.append(dict(T_reheat_GeV=x["T_reheat_GeV"],omega_scalar_from_central_parent_source=s,
               unavoidable_omega_stable_from_the_same_parents_if_conserved=stable_parent,
               separate_direct_lepton_to_stable_pair_omega_central=stable_light,
               fractional_stable_loss_required_to_match_total_with_scalar_held=1-fraction,
               if_every_stable_converts_one_for_one_scalar_omega_central=2*s+y,
               if_every_stable_converts_with_parent_envelope_and_central_light_scalar_omega=2*upper+y))
    return dict(scope="Matched prompt daughter identity and conditional stable-only two-to-two number ledger. "
                      "Direct light source is central vacuum charged-lepton production, integrated down to 1 MeV; "
                      "no pion-to-light, heavy regeneration, independent production, late entropy, or parent re-entry. "
                      "All-conversion comparison grants one scalar per stable without a rate or a momentum model. "
                      "Parent enclosure plus central light is not an upper for arbitrary light production.",
                m_stable_over_scalar=ratio,omega_dark_reference=total,
                prompt_only_scalar_mass_share_if_counts_conserved_upper=1/(1+ratio),rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos",type=Path)
    p.add_argument("bracket",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.bracket),indent=2,allow_nan=False))
