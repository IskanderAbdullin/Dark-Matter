# A live, number- and energy-conserving scalar exterior in a decay-only control

Start both arms from the same final hard and cold state of the [accepted 216-step fine trajectory](INDEPENDENT-EARLY-POLE-FINE-EXTENDED-MEMORY-2026-11-27.md). In the second arm, move positive archived spectator counts into live scalar variables at their existing native momentum nodes. Keep the original coupling, scalar gap equation, hard moments and initial total number and energy. The promoted volume uses the native shell volume minus any overlap already represented by the cold window. This changes which modes evolve, not the initial distribution. Route each collision through the [reciprocal outer packet ledger](INDEPENDENT-EARLY-POLE-SCALAR-EXTERIOR-RECIPROCAL-LEDGER-2026-11-27.md). The short continuation includes only decay and inverse decay, not the other hard channels from the original trajectory or expansion.

| Promoted node support and hard reconstruction | Extension exposure | Integrated exterior gain |
| --- | ---: | ---: |
| q=.1–5, native steps (11 modes) | 2e−6 | 5.496615e−13 |
| q=.1–5, count-conserving continuous (11 modes) | 2e−6 | 5.522148e−13 |
| q=.1–5, native steps | 7.2e−5 | 1.97874e−11 |
| q=.1–5, count-conserving continuous | 7.2e−5 | 1.98793e−11 |
| q>1e−5 through 5, native steps (37 modes) | 2e−6 | 8.779168e−13 |

The initial cold rate agrees across arms (the native value is .1177324264); the initial scalar gap agrees within 1.1e−19 and the particle partition agrees to about 4e−18. All cases retain the signed parent/daughter count, energy and entropy ledgers. The coarse native one-node warm-shell representation gives an initial q≥.1 source about 5% larger than the preceding continuous momentum fold. This is a momentum-discretization limitation, not an enhancement from feedback.

A separate centered initial-direction test isolates the effect of adding the exterior packets on the instantaneous cold source. With only the q=.1–5 nodes the added derivative is about −3.9e−7 in either hard reconstruction. Allowing the whole high-q exterior in the native-step state makes it about −1.50e−5; at the largest diagnostic offset the hard-occupation direction accounts for −1.44e−5 and the scalar pole shift for about −6.31e−7. At short exposure 2e−6, the latter combined derivative would give a local quadratic correction of order −3e−17; the independently integrated broad-exterior runs give differences of that order. The q=.1–5 warm gain at the longer exposure is stable to better than 1e−11 relative when the maximum time step is halved. The much smaller cold-gain difference, however, changes sign or amplitude between these controls, so separate adaptive trajectories do not resolve it quantitatively. A converged late abundance or an observational hot-component fraction cannot be inferred.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
b=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
for p in steps cell-conservative; do
 python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" \
  --projection "$p" --duration 2e-6 --panels 8 --radial-order 12
 python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" \
  --projection "$p" --duration 7.2e-5 --panels 8 --radial-order 12
 # Repeat the preceding 7.2e-5 control with --panels 16.
 python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" \
  --projection "$p" --linear-probe --radial-order 12
done
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --outer 1.00001e-5 5
python -W error early_pole_live_scalar_exterior_decay_reference.py "$b" --outer 1.00001e-5 5 \
 --panels 8 --radial-order 12 --linear-probe
```
Raw: `results/early-pole-live-scalar-exterior-decay-2026-11-27.json.gz`.
