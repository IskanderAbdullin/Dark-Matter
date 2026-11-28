"""Opt-in CLASS output-only omission of very large exact-mode history tables."""
import argparse
import hashlib
import json
from pathlib import Path

BASE='if (pop->write_perturbations == _TRUE_ && ppt->has_perturbations) {'
UPDATED=BASE.replace('has_perturbations)', 'has_perturbations && getenv("CLASS_EXACTK_POWER_ONLY") == NULL)')


def run(source,patched):
    assert source.resolve() != patched.resolve()
    data=source.read_bytes()
    assert data.count(BASE.encode())==1 and b'CLASS_EXACTK_POWER_ONLY' not in data
    out=data.replace(BASE.encode(),UPDATED.encode())
    patched.write_bytes(out)
    dig=lambda d:hashlib.sha256(d).hexdigest()
    return dict(scope="After the solver has populated output arrays, optionally skip writing perturbation tables. "
                      "This does not bypass integrating/storing exact modes or change memory needs for those arrays.",
                original_source_sha256=dig(data),patched_source_sha256=dig(out))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source",type=Path)
    p.add_argument("patched",type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.source,a.patched),indent=2))
