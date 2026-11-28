"""Certify the selected direct-parent plus direct-light number budget under restricted reactions."""
import argparse
import hashlib
import json
from pathlib import Path

from archived_collisionless_streaming_reference import Background
from three_species_decay_onsager_reference import M, m, mp
from warm_direct_vacuum_absolute_source_envelope_reference import bounds, log_slope_minimum, read


def run(eos,kloe,gap,ladder):
    b=Background(eos)
    a=read(ladder)
    inputs=read(gap)
    scale=a["Y_to_Omega_h2"]
    slopes=[log_slope_minimum(b.gs,3),log_slope_minimum(b.gr,4)]
    out=[]
    for row in a["rows"]:
        top=row["T_reheat_GeV"]
        results=[]
        for n in (240,480):
            parent=bounds(b,kloe,inputs,n,top,mass=M)
            light=bounds(b,kloe,inputs,n,top,mass=m)
            p=parent["total_known_plus_unit_high_yield_upper"]*scale
            l=light["total_known_plus_unit_high_yield_upper"]*scale
            prev=next(v["upper"] for v in row["controls"] if v["invariant_panels"]==n)
            assert abs(p/prev-1.)<1e-13
            V=2*p+l
            results.append(dict(panels=n,heavy_pair_scalar_equivalent_upper=p,
                                direct_light_scalar_equivalent_upper=l,
                                total_count_scalar_equivalent_upper=V,
                                fraction_of_reference_half_dark_upper=V/.060,
                                full_parent_source_bounds=parent,full_light_source_bounds=light))
        assert results[1]["total_count_scalar_equivalent_upper"]<results[0]["total_count_scalar_equivalent_upper"]
        assert results[1]["heavy_pair_scalar_equivalent_upper"]<results[0]["heavy_pair_scalar_equivalent_upper"]
        assert results[1]["direct_light_scalar_equivalent_upper"]<results[0]["direct_light_scalar_equivalent_upper"]
        out.append(dict(T_reheat_GeV=top,controls=results))
    assert out[-1]["controls"][-1]["total_count_scalar_equivalent_upper"]<.060
    assert all(x["controls"][-1]["total_count_scalar_equivalent_upper"]<
               y["controls"][-1]["total_count_scalar_equivalent_upper"] for x,y in zip(out,out[1:]))
    return dict(scope="Initially empty, prescribed adiabatic visible radiation through .1 MeV. "
                      "Direct selected vacuum portals to heavy and light equal-mass pairs have the same supplied "
                      "vector coupling. Charged leptons and the charged-pion measured window, two lower-band "
                      "positive caps and assumed R_pi<=1 in the high band are used in both. "
                      "Individual heavy decay/reverse preserves the weighted count Q=N_scalar+N_light+2N_heavy. "
                      "Stable pair conversion/reverse preserves Q; forward heavy-pair to scalars only lowers it. "
                      "Do not use this if inverse heavy pair creation from dark light particles or other "
                      "unrepresented number sources may increase Q. No depletion rates or transfer predicted.",
                masses_GeV=dict(parent=M,stable=m,scalar=mp),
                eos_sha256=hashlib.sha256(eos.read_bytes()).hexdigest(),
                KLOE_sha256=hashlib.sha256(kloe.read_bytes()).hexdigest(),
                eos_min_slopes=slopes,rows=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("eos","kloe","gap","ladder"):
        p.add_argument(key,type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kloe,a.gap,a.ladder),indent=2,allow_nan=False))
