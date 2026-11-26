# A conservative reverse-coordinate partition removes an isolated cold-cell spike

The [tagged finite-momentum audit](INDEPENDENT-WARM-COLD-LIGHT-FINITE-PAIR-2026-11-26.md)
found that projecting scalar-pair reverse events from samples of both initial
scalars severely undersampled light particles in the first cell near 8.133 MeV.
The same native projection instead **overcounts** that inverse gain on the
7.786–7.755 MeV frozen snapshots, then almost misses it again around 7.724 MeV.
This is not a stable multiplicative correction to a genuine physical channel.

On each inverse event, partition by whether at least one outgoing fermion
falls into the first light cell. Keep the native incoming-scalar coordinates
for the complement. Integrate the cold sector in tagged final-light variables
and scatter the same event to both scalar losses and both fermion gains.
Because either outgoing charge can tag an event, divide its weight by two
if both light particles are in the cell. The forward operator is unchanged.
Projecting all four legs conserves count and energy per event; the frozen
pair and the evolved block both check their ledgers independently.

On a synthetic broad-cell control, about half of the restricted events have
both light products in the cell. There the two independent coordinate
choices agree on the four projected distributions after cut subtraction
to 8.5e-4 in normalized shell L1 at the finer paired rule. In the actual q64 bath, only about
5.4–5.8e-9 of inverse events hit the first cell under the hybrid rule, and
the two-cold subset is far smaller. A separate scan of 12 frozen q64 shells
shows the local native quadrature error extends beyond the first cell;
the current replay test repairs only the first one.

Replay the accepted q32 block to 7.754805 MeV from the same two hashed
predecessors, with the same FRW step, operators and time grid. The native
result reproduces the archived endpoint hash. The native cell mean rises
at +1.43547e-4 per unit ln(a); with the conservative partition it falls at
−4.71847e-5, with the same sign and close agreement on both tagged
quadrature levels. The changes to the complete scalar and light shell
distributions are below 3.4e-10 in all three normalized moment L1 tests
(number, energy and pressure) for this same block. This establishes a
local numerical spike, not a correction to the full spectrum history,
a convergence proof for adjacent cold cells, or a relic prediction.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p replay/q32 replay/q64
for q in 32 64; do
  for j in 10 11 12; do
    zcat "results/warm-cold-pair-spike-q$q-step$j-2026-11-26.json.gz" > "replay/q$q/step-$j.json"
  done
done
python -W error warm_cold_pair_spike_replay_reference.py EOS.dat replay replay/32-low --q 32
python -W error warm_cold_pair_spike_replay_reference.py EOS.dat replay replay/32-fine --q 32 --fine
python -W error warm_cold_pair_partition_controls_reference.py \
  results/warm-recoil-crossing-q64-terminal-2026-11-26.json.gz \
  results/warm-cold-pair-spike-q64-step12-2026-11-26.json.gz
python -W error warm_cold_light_finite_pair_reference.py EOS.dat --snapshots \
  results/warm-cold-pair-spike-q64-step{11,12,13}-2026-11-26.json.gz
~~~

[Matched q32 replays](results/warm-cold-pair-spike-replays-2026-11-26.json);
[full coordinate and shell controls](results/warm-cold-pair-partition-controls-2026-11-26.json);
[reported uncorrected q64 intermediate-rate checks](results/warm-cold-pair-spike-q64-frozen-2026-11-26.json).

Exact [refined lifted q32 output](results/warm-cold-pair-spike-q32-lifted-2026-11-26.json.gz).

## Separate time-step control

Same q32 predecessors, endpoint temperature, quadrature and partition; subdivide
only the original expansion interval. Both schemes keep opposite signs across
the tested time ladder. Rates are cell-mean occupation per unit ln(a).

| matched cuts | original | lifted |
|---:|---:|---:|
| 1 | +1.43547e-4 | -4.71847e-5 |
| 2 | +1.257074e-04 | -4.771821e-05 |
| 4 | +1.209726e-04 | -4.799308e-05 |

~~~sh
for n in 2 4; do
  python -W error warm_cold_pair_spike_replay_reference.py EOS.dat replay replay/32-time-$n --q 32 --subdivisions "$n"
done
~~~

[Time-step diagnostics](results/warm-cold-pair-spike-time-2026-11-26.json).
