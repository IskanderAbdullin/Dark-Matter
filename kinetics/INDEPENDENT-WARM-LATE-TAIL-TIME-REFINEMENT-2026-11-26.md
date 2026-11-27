# Follow the cut-aligned .68 GeV excess through the first reversal and decay

Keep the exact-vacuum initial state and the .68 GeV
choice of the [cut-aligned prompt test](INDEPENDENT-WARM-CUT-ALIGNED-PROMPT-VACUUM-CONTROL-2026-11-26.md).
The branch uses a fixed clock, a closed 1.05 GeV momentum
pad and only the prompt operator. Compare the original
three accepted implicit microsteps per reported interval
against six per interval **after the shared, exact three-step
prefix to Hdt=1e-11**. The first interval is not a refinement
test. Starting the six-step ladder there avoids substituting
an artificial vacuum floor or loosening the entropy guard.

Both complete time ladders to Hdt=1e-6 preserve the
positive-to-negative sign change between 1e-8 and 1e-7.
The maximum difference in sampled excess over the larger
sampled excess is 3.74%. The small late residual is more
sensitive in relative terms. Define total energy above the
cut as the full energy of every above-cut daughter. It
is distinct from the excess, which subtracts the cut
energy for each such daughter.

| Hdt | Three-step excess / initial energy | Shared-prefix six-step excess / initial energy | Six-step total energy above cut / initial energy |
| ---: | ---: | ---: | ---: |
| 1e-8 | 1.102e-18 | 1.130e-18 | 1.211e-16 |
| 1e-7 | 2.275e-18 | 2.326e-18 | 2.447e-16 |
| 1e-6 | 2.701e-19 | 1.834e-19 | 1.870e-17 |

The excess shrinks by at least 88% between the last two
samples in both ladders but does not vanish. This is a
resolved state-moment result for one chosen initial shape,
not a global limit on the lifetime or abundance of a tail.

## Lock and inspect the actual steps near the reversal

From the same accepted six-step state at Hdt=1e-8,
replay only the interval to 1e-7 using locked uniform
six- and twelve-step clocks. Accept no silent subdivision.
At each accepted root verify the exact discrete identity
`W_new−W_old = h*F_new` for the cut-aligned excess moment,
where h is one Hubble-normalized microstep. The signed
sum closes within 4e-32 in units of the initial dark energy.

The six-step sampled sign bracket is [7.0e-8, 8.5e-8];
the twelve-step bracket is [6.25e-8, 7.0e-8]. Narrower
sampling and smaller implicit time steps shift the bracket.
These are intervals with opposite signs on the two
numerical trajectories, not an enclosure of the continuum
turning time. The independent angular emission term falls
from 1.415e-10 to around 9.7e-11 in retained-energy Hubble
units, while angular consumption grows from 5.128e-11 to
around 1.03e-10. The reversal is the change of balance of
these two terms, not the gross inverse source shutting off.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/warm-corrected-anchor-signed-floor-2026-11-26.json
b=results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz
s=results/warm-late-cut-aligned-hard068-test-step1-2026-11-26.json.gz
python -W error warm_corrected_late_prompt_cut_aligned_reference.py EOS.dat "$a" "$b" n3 --hard .68 --substeps 3 --predictor
python -W error warm_corrected_late_prompt_cut_aligned_compact_after_prefix_reference.py EOS.dat "$a" "$b" n6 --hard .68 --substeps 6 --predictor --prefix "$s"
python -W error warm_late_cut_aligned_prompt_clock_audit_reference.py n3 n6
python -W error warm_late_cut_absolute_tail_energy_reference.py "$a" n3 n6
for n in 6 12; do
  python -W error warm_late_cut_turning_clock_reference.py EOS.dat "$a" n6/step-4.json "turn$n" --substeps "$n" >"turn$n.json"
done
python -W error warm_late_cut_turning_convergence_reference.py turn6.json turn12.json
~~~

[Full time comparison](results/warm-late-aligned-hard068-time-refinement-2026-11-26.json.gz);
[absolute and excess energies](results/warm-late-aligned-hard068-total-tail-energy-2026-11-26.json.gz);
[locked turning diagnostics](results/warm-late-aligned-turning-six-twelve-2026-11-26.json);
[raw six-step turning tape](results/warm-late-aligned-turning-six-2026-11-26.json.gz);
[raw twelve-step turning tape](results/warm-late-aligned-turning-twelve-2026-11-26.json.gz).

[Terminal three-step archive](results/warm-late-aligned-hard068-terminal-3-2026-11-26.json.gz);
[terminal shared-prefix six-step archive](results/warm-late-aligned-hard068-terminal-6-2026-11-26.json.gz).
[Exact shared turning seed](results/warm-late-aligned-turning-seed-2026-11-26.json.gz).
