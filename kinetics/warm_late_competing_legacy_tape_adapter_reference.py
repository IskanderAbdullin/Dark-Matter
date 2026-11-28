"""Bind the same fixed physical cut to the older cut-aligned late tape."""

import argparse
from pathlib import Path

import numpy as np

from three_species_decay_onsager_reference import M,m,mp
from warm_coupled_serial_restart_reference import atomic,sha
from warm_light_occupation_stage_audit_reference import digest,read


def adapt(prior,source,out):
    dossier=read(prior)
    cut=float(np.hypot(M,1./dossier["branches"][-1]["entropy_expansion_from_original_anchor"]))
    out.mkdir(exist_ok=True,parents=True)
    for j in (0,4,5,6):
        f=source/f"step-{j}.json"
        z=read(f)
        S=z["snapshot"]
        assert z["hash"]==digest(S) and z["plan"]["artificial_empty_shell_filling"]==0.
        target=[np.sqrt((cut-d)*(cut+d)) for d in (mp,m)]
        q=np.asarray(S["lab_momentum_edges_GeV"])
        assert all(np.any(q==edge) for edge in target)
        a=dict(z,plan=dict(z["plan"],normalized_cut_GeV=cut,physical_cut_source_sha256=sha(prior)),
               signed_native_flux=z["native_signed_cut_flux"],raw_tape_file_sha256=sha(f))
        atomic(out/f.name,a)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("prior","source","output"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    adapt(a.prior,a.source,a.output)
