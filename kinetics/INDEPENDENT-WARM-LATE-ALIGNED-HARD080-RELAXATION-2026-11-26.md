# Follow the initially above-cut .80 GeV cohort realization

Repeat the [cut-aligned .68 experiment](INDEPENDENT-WARM-LATE-TAIL-TIME-REFINEMENT-2026-11-26.md)
with the .80 GeV hard packet. All 536 receipts still preserve
their number and birth-energy mean, but this realization
already has substantial *relative* excess above the energy
cut. The cut, momentum pad (1.05 GeV), selected closed
prompt channel and radiation clock remain prescribed.
No vacuum occupancy floor is inserted.

Both time ladders share the exact three-microstep prefix
to Hdt=1e-11; thereafter compare three and six microsteps
per reporting interval. The signed excess flux is negative
at **all saved points** on both ladders, unlike the .68 path.
No interval between samples is certified sign-definite.
The earlier upper bounds on the instantaneous positive source
are not used as a lifetime estimate.

| Observable, normalized to initial represented dark energy | Start | Three steps per later interval at Hdt=1e-6 | Six steps per later interval at Hdt=1e-6 |
| --- | ---: | ---: | ---: |
| Convex excess above the cut | 9.845e-14 | 2.034e-14 | 1.814e-14 |
| Full daughter energy above cut | 9.240e-13 | 2.006e-13 | 1.796e-13 |

The two runs remove approximately 79% and 82% of the initial
excess over the reported window. The maximum difference in
excess at any common sample is 2.24% of the initial value;
the terminal residual differs by about 12% relative to the
finer one. The terminal current remains negative
(−2.884e-8 per Hubble time in units of the initial energy
for the finer run). Neither the late value nor the slope
is fitted as a universal exponential or extrapolated beyond
the last accepted state.

This comparison explains how two matched-moment initial shapes
can produce different amounts and directions of early transfer.
It does not evolve visible, other dark or expansion effects,
or the unknown physical exterior distribution.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
s=results/warm-late-aligned-hard080-prefix-2026-11-26.json.gz
python -W error warm_corrected_late_prompt_cut_aligned_reference.py EOS.dat "$a" "$b" n3 --hard .80 --substeps 3 --predictor
python -W error warm_corrected_late_prompt_cut_aligned_compact_after_prefix_reference.py EOS.dat "$a" "$b" n6 --hard .80 --substeps 6 --predictor --prefix "$s"
python -W error warm_late_cut_aligned_prompt_clock_audit_reference.py n3 n6
python -W error warm_late_cut_absolute_tail_energy_reference.py "$a" n3 n6
~~~

[Exact shared prefix](results/warm-late-aligned-hard080-prefix-2026-11-26.json.gz);
[complete paired trajectories](results/warm-late-aligned-hard080-time-refinement-2026-11-26.json.gz);
[above-cut energy and excess](results/warm-late-aligned-hard080-total-tail-energy-2026-11-26.json.gz);
[coarse terminal](results/warm-late-aligned-hard080-terminal-3-2026-11-26.json.gz);
[fine terminal](results/warm-late-aligned-hard080-terminal-6-2026-11-26.json.gz).
