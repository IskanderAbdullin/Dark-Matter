# Table-refine the moderate-occupation path and locate the end of fast heat sharing

Repeat the [early path to 40 MeV](INDEPENDENT-WARM-REHEATING-MODERATE-OCCUPATION-2026-11-26.md)
with FD/BE thermodynamics and iterated numerical quantum
rate tables. Each iteration recalculates rate factors on
the updated path. The light–scalar branch uses direct
initial and final occupation integrals; the thermal
reverse rate uses the exact chemical-affinity relation.
The heavy forward branches retain the midpoint of the
local final heavy Pauli interval as an explicit approximation.

| Early rate table | Sampled max scalar phase factor | Scalar invariant at 40 MeV |
|---:|---:|---:|
| 11 path nodes / coarser angular grid | .638094 | 2.21162e−4 |
| 16 path nodes / finer angular grid | .638329 | 2.21204e−4 |

The .65 test ceiling is not reached in either arm.
The maximum is sampled near 62.44 MeV, not an exact
continuous maximum. The independent higher-occupation
light–scalar rate at the baseline peak is 1.02799 times
its matched-grid MB value; at 40 MeV it is 1.01850 times.
Take the FD/BE and MB ratio on the same positive nodes to
avoid confusing a small occupation correction with the
difference between independent quadratures. Both measure
the same sampled thermal states and neither covers an
arbitrary nonthermal distribution.

## Cooler diagnostic from the MB-event-rate handoff at 40 MeV

Start this separate control from the MB-rate arm, retain
FD/BE thermodynamics and revert to tabulated MB event
rates. The sampled saturated normal-branch capacity
fraction is about .838 at 18.37 MeV and .736 at 8 MeV.
At 8 MeV the thermal fast heat gap remains about
10,359 times H, or about 6,397 times the larger of H
and the selected fractional number rates. When continued
with a declared clock-ratio stop of 30, both numerical
grids stop at visible 6.000 MeV. The dark temperature
there is 2.485 MeV and the light/scalar invariants are
5.021e−5 and 6.091e−4, respectively. This is a switch
criterion for a model in which the two species have shared
one temperature, not a physical temperature at which
a process is known to end. Heavy decays and inverse rates
cannot be inferred to remain fast indefinitely.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_tabulated_quantum_reference.py EOS.late.dat --early --wide
python -W error warm_prompt_reheating_early_seed_sensitivity_reference.py EOS.late.dat --peak .65 > mids.json
python -W error warm_prompt_reheating_thermal_late_probe_reference.py EOS.late.dat mids.json > late.json
python -W error warm_prompt_reheating_thermal_late_probe_reference.py EOS.late.dat mids.json --cold > cold.json
python -W error warm_prompt_reheating_quantum_heat_gap_reference.py EOS.late.dat late.json
~~~

[Quantum table summary](results/warm-prompt-reheating-wide-quantum-summary-2026-11-26.json);
[direct moderate-occupation rates](results/warm-prompt-reheating-mid-thermal-rates-2026-11-26.json);
[direct low-temperature rate controls](results/warm-prompt-reheating-deep-thermal-rates-2026-11-26.json);
[8 MeV path](results/warm-prompt-reheating-late-probe-2026-11-26.json.gz),
[heat gap](results/warm-prompt-reheating-deep-heat-2026-11-26.json.gz),
[clock-limited extension](results/warm-prompt-reheating-cold-probe-2026-11-26.json.gz).
