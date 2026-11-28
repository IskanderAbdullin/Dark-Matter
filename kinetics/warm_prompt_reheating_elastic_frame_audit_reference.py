"""Independent trace, detailed-balance, and cutoff tests of the elastic heat fold."""
import argparse
import json
from pathlib import Path
import numpy as np

from scalar_rest_yukawa_elastic_ir_reference import cross
from three_species_decay_onsager_reference import m,mp
from warm_prompt_reheating_two_temperature_contact_probe_reference import exact_elastic_sigma_at_angle,fold


def run(path):
    S=np.geomspace(.000401,.003,20)
    z=np.linspace(-1.,1.,56)
    rest=np.array([cross(1.,(x-m*m-mp*mp)/(2*mp),-z) for x in S])
    inv=exact_elastic_sigma_at_angle(S[:,None],z)
    err=np.max(abs(rest/inv-1))
    assert err<5e-12
    d=json.loads(path.read_text())
    b=next(a for a in d["controls"][-1]["heat_arms"] if a["prompt_heat"])["rows"]
    Q,M,X=fold(.002,.002,[.014,.3],(8,12,13),True)
    ratio=abs(Q)/(M*.002)
    assert ratio<1e-14
    comparisons=[]
    for i in (4,16):
        r=b[i]
        Tl,Ts=r["light_and_scalar_dark_T"]
        pk=r["light_and_scalar_phase_peak"]
        result=[]
        for cut,az in ((.0024,10),(.003,10),(.003,19)):
            q,e,out=fold(Tl,Ts,pk,(9,13,16),True,cut,az)
            result.append(dict(s_max=cut,azimuth_count=az,off_shell_heat_GeV5=q,
                               gross_below_cut_event_density_GeV4=e,
                               outside_contact_proxy_share=out/(out+e)))
        assert abs(result[-2]["off_shell_heat_GeV5"]/result[-1]["off_shell_heat_GeV5"]-1)<2e-6
        comparisons.append(dict(T_visible_GeV=r["visible_T"],controls=result))
    return dict(scope="Selected below-pole invariant angle only; both cuts strictly exclude the pole. "
                      "Contact-only proxy outside each cut is not a bound on omitted near-pole physics.",
                invariant_vs_rest_trace_relative_defect=float(err),
                common_temperature_zero_heat_over_gross_x_T=ratio,case_checks=comparisons)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("path",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.path),indent=2,allow_nan=False))
