# The dilute continuation depends on when it is started

The [signed model from the 60 MeV seed](INDEPENDENT-WARM-REHEATING-EXPANDING-SLOW-2026-11-26.md)
starts from an imposed distribution with no earlier
dark reaction history. Test a distinct, early-seeded
continuation for the same 100 MeV instantaneous
reheating. Initialize each arm at a different bath
temperature just below reheating with the direct
lepton-only seed moments available at that time.
Then evolve the same selected signed MB collision
rates with FD/BE thermodynamic moments, expansion
and continuing visible lepton input.

| Initial seed (MeV) | Visible T at phase factor .075 (MeV) | Visible T at phase factor .28 (MeV) | Scalar invariant yield at the .28 stop |
|---:|---:|---:|---:|
| 99.98 | 90.7514 | 86.0943 | 7.06060e−5 |
| 99.90 | 90.7514 | 86.0942 | 7.06060e−5 |
| 99.50 | 90.7503 | 86.0938 | 7.06074e−5 |
| 99.00 | 90.7429 | 86.0909 | 7.06169e−5 |
| 97.00 | 90.5520 | 86.0147 | 7.08682e−5 |

The early arms agree closely on their respective
stop temperatures, but do **not** stay dilute to
60 MeV in this reduced model. The preset 60 MeV seed
cannot be described as the endpoint of these early
controlled trajectories. The phase factors are
numerical applicability markers, not measured
onsets, condensation thresholds, or late relics.

The common-temperature assumption also needs a
check close to vacuum. Match thermal numbers and
energy to the direct FD input moments at temperatures
only a few keV below reheating. Compute the full FD/BE
linear heat conductance of the fast prompt operator
using the final-state weighted variance at fixed parent
energy, and compare it to the source rate divided by
the inherited count. The ratio is .927 at 2 keV below
reheating and 1.448 at 2.5 keV below; at 20 keV below,
the ratio is 92.7. At the earliest followed seed, the
heat exchange is already faster than the source by this
local test. This is not a proof that all momentum modes
have equilibrated in the initial boundary interval.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_early_seed_sensitivity_reference.py EOS.late.dat > seeds.json
python -W error warm_prompt_reheating_onset_heat_reference.py EOS.late.dat
python -W error warm_prompt_reheating_quantum_heat_gap_reference.py EOS.late.dat seeds.json
~~~

[Seed ladder](results/warm-prompt-reheating-early-seed-sensitivity-2026-11-26.json);
[initial heat/source race](results/warm-prompt-reheating-onset-heat-2026-11-26.json);
[heat-rate and finite-difference cross along the earliest path](results/warm-prompt-reheating-early-path-heat-2026-11-26.json).

A [separate direct quantum rate refinement](INDEPENDENT-WARM-REHEATING-EARLY-TABULATED-QUANTUM-2026-11-26.md)
confirms the earlier applicability stop with the same
early seed, instead of continuing directly to 60 MeV.
