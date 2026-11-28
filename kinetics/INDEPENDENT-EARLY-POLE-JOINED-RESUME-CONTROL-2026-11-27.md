# Checked joined-model continuation in its native solver coordinates

Long charged-logit simulations can resume a finished segment
without reconstructing the charged cell occupations from reported
moments. The optional joined solver switches
`--emit-checkpoint` and `--checkpoint-input` save and reload the
full native chart: hard coordinates, cold occupations, accumulated
ledger variables, and the absolute exposure at the end of a segment.
They are limited to the native and exact archive-replacement arms.

The resumed program checks its state schema and the recorded model
provenance before solving. The manifest records SHA-256 hashes of the
initializer, top-level program, imported model Python sources and
selected spectra, plus Python/NumPy/SciPy and long-double precision.
Projection, archive intervention, hard depth, elastic choice,
tolerances and the maximum nominal step must match. Reconstructed
global energy, number and entropy/ledger constraints and positive
cold occupations are rechecked before resuming. This is a numerical
compatibility gate, not a signature authenticating a hand-edited
state.

The retained **H1 local-daughter cold32** control compares an
uninterrupted total exposure 2.666666666666667e−9 with two consecutive
segments of half that interval. The resumed cold join is exactly
the saved boundary array. The largest terminal mode difference
is 4.34e−19 and the largest terminal native chart difference
is 9.78e−13 after the integrator's step history restarts.
The optional uninterrupted run reproduces every ordinary report field
of the prior no-checkpoint solver exactly on this matched fixture.
An altered program fingerprint and a negative cold state are both
rejected. These are short-regime tests, not validation for long H8
histories or a promise of step-for-step identity.

Example from kinetics; write a partial response to a temporary file
before making it a checkpoint:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
seed=results/early-pole-first-window-hard1-cold32-paired-initializer-2026-11-27.json.gz
script=early_pole_joined_simultaneous_chart_reference.py
python -W error "$script" "$seed" --local-from-paired-archive --duration 1.3333333333333335e-9 --panels 1 --record-mode-history --emit-checkpoint > first.json.tmp
mv first.json.tmp first.json
python -W error "$script" "$seed" --local-from-paired-archive --duration 1.3333333333333335e-9 --panels 1 --record-mode-history --emit-checkpoint --checkpoint-input first.json
```

Quickly check the preserved three-way receipt:
```sh
python -W error early_pole_joined_resume_receipt_reference.py \
  results/early-pole-joined-resume-local-cold32-raw-control-2026-11-27.json.gz
```
[Summary](results/early-pole-joined-resume-local-cold32-2026-11-27.json).
