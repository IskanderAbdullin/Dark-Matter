"""Checkpoint a warm continuation while lifting inverse events over a cold cut."""

import argparse
import json
from pathlib import Path

import coupled_warm_pair_continuation_reference as kernel
import warm_cold_light_pair_partition_reference as pair
from warm_coupled_serial_restart_reference import atomic, sha
from warm_light_occupation_stage_audit_reference import read
from warm_q64_history_serial_reference import run


def execute(eos, tape, folder, step, count, cut, tag_order):
    folder.mkdir(parents=True, exist_ok=True)
    # The inner warm runner already checks its original kernel, input and
    # clock. This extra manifest pins the substituted collision and options.
    manifest = folder / "cold_partition_manifest.json"
    local = Path(__file__).resolve().parent
    modules = [Path(__file__), Path(pair.__file__),
               local / "archived_scalar_pair_reaction_reference.py",
               local / "warm_cold_light_finite_pair_reference.py"]
    keys = dict(cut_cells=cut, tagged_order=list(tag_order), relative_step=step,
                steps=count, seed_sha256=sha(tape), eos_sha256=sha(eos),
                module_sha256={m.name:sha(m) for m in modules})
    if manifest.exists():
        assert read(manifest) == keys
    else:
        assert not list(folder.glob("step-*.json"))
        atomic(manifest, keys)
    saved = kernel.rate
    peak = 0.

    def rate(p, th, order):
        nonlocal peak
        V, err, ratio, diagnostics = pair.split_rate(p, th, order, tag_order, cut_cells=cut)
        peak = max(peak, diagnostics["tag_closure"])
        return V, err, ratio

    try:
        kernel.rate = rate
        R = run(eos, tape, step, count, folder)
    finally:
        kernel.rate = saved
    return dict(scope="Selected warm kernel with a conservative cold reverse partition, no unselected channels.",
                manifest=keys, max_copy_symmetry_defect_this_execution=peak, continuation=R)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("eos", "tape", "folder"):
        p.add_argument(key, type=Path)
    p.add_argument("--step", type=float, required=True)
    p.add_argument("--blocks", type=int, required=True)
    p.add_argument("--cells", type=int, choices=(1,3,6,9), default=9)
    p.add_argument("--fine", action="store_true")
    a = p.parse_args()
    order = (6,12,9,12,12) if a.fine else (4,8,6,8,8)
    print(json.dumps(execute(a.eos,a.tape,a.folder,a.step,a.blocks,a.cells,order),
                     indent=2, allow_nan=False))
