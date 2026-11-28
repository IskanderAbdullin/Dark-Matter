# Extend the matched warm states and inspect actual light occupations

Resume the two complete [32-block q16 and q32 checkpoints](INDEPENDENT-WARM-EARLY-SUPPORT-32-BLOCKS-2026-11-26.md)
at 10.442 MeV, preserving the positive exterior restarts.
Evolve on four new logarithmic expansion blocks from that
anchor through b=1.02, with T just below 10.238 MeV. This is
not the separate longer q32 outcome. Each stage is stored
and rechecked against the parent snapshot, plan hash and
accumulated ledger. No internal steps were rejected.

At each saved state independently evaluate the largest
within-shell light and scalar occupations from the linear
entropy endpoints. Also integrate the stationary heavy
energy for decay/inverse alone with the daughters frozen,
using two radial orders. This is not a new coupled heavy
evolution.

| q | Initial → final retained light peak | Initial → final retained scalar peak | Initial → final frozen heavy energy fraction |
| --- | ---: | ---: | ---: |
| 16 | 6.07630×10⁻⁵ → 5.95463×10⁻⁵ | .000375252 → .000380333 | 1.42517×10⁻⁶ → 1.23339×10⁻⁶ |
| 32 | 5.86644×10⁻⁵ → 5.79567×10⁻⁵ | .000376179 → .000381167 | 1.42900×10⁻⁶ → 1.23688×10⁻⁶ |

The light peak decreases at each saved time. Tracked
retained-plus-new-exterior comoving number is conserved
to (2×10^{-14}) relative across these saved states.
This does not bound intervening times, omitted earlier
exterior occupancies or the full warm evolution.
These finite data cannot validate the
[analytic per-spin cap](INDEPENDENT-WARM-SIGNED-FERMION-TIME-COMPARISON-2026-11-26.md)
for all times and all momenta.

The tapes contain each of the four full atomic results.
Construct the source as the previous snapshot wrapped
in a one-entry snapshots object using the driver's
bytes_of serializer. Invoke the same warm serial driver
with --stop 1.02 --steps 4 --local-prompt-jacobian
--max-blocks 1 on each of four passes, saving each
atomic checkpoint and retaining the same state directory.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_light_occupation_stage_audit_reference.py \
 "$D"/warm-early-q{16,32}-time32-checkpoint-2026-11-26.json.gz \
 "$D"/warm-light-occupation-next-q{16,32}-tape-2026-11-26.json.gz
~~~
