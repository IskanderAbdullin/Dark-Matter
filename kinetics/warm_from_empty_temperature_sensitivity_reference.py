"""Map where the restricted vacuum upper ceases to exclude half of 0.120."""
import argparse
import hashlib
import json
from pathlib import Path

from warm_from_empty_recycling_envelope_reference import with_eos


def run(eos, kernel):
    K = json.loads(Path(kernel).read_bytes())["simultaneous_energy_density_coefficient_GeVminus2"]
    controls = []
    previous = -1.
    for TR in (.0545,.055,.0555,.056):
        low, high = [with_eos(TR,K,n,d,eos) for n,d in ((700,.01),(1300,.006))]
        y = high["Omega_phi_total_ceiling"]
        assert y < low["Omega_phi_total_ceiling"] and y > previous
        previous = y
        controls.append(dict(T_R_GeV=TR,coarse=low,refined=high,
                             fine_upper_excludes_half=bool(y<.060)))
    return dict(half_reference=.060,controls=controls,sha256=[
                    hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in (eos,kernel)],
                scope="Only the numerical upper on vacuum charged-lepton instantaneous "
                      "production. Crossing the reference is loss of exclusion by the "
                      "upper, not a valid solution or a physical threshold.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("eos")
    p.add_argument("kernel")
    a=p.parse_args()
    print(json.dumps(run(a.eos,a.kernel),indent=2,allow_nan=False))
