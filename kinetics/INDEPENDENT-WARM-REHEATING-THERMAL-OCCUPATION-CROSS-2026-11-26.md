# Directly integrate occupation factors on three states of the expanding comparison

The [three expanding sensitivity trials](INDEPENDENT-WARM-REHEATING-EXPANDING-SLOW-2026-11-26.md)
use local endpoint bounds on final and initial occupancy
factors, not integrated quantum rates. At three states of
the midpoint trial, hold T and chemical potentials fixed
and evaluate the light pair ↔ scalar pair integrals with
full FD/BE factors. Compare the directly evaluated reverse
and forward event ratio against the exact thermal affinity
and repeat with a finer momentum and angle sampling.

| Dark temperature (MeV) | Light-to-scalar forward quantum multiplier over MB | Relative defect in reverse/forward balance |
|---:|---:|---:|
| 54.3269 | 1.000888 | 1.91e−6 |
| 34.4589 | 1.001302 | 6.55e−7 |
| 24.3404 | 1.002456 | 3.27e−7 |

At the latest tested state, the available positive bracket
was approximately [0.939,1.169], much wider than the actual
occupation correction for this cross section on this state.
The independent light- and scalar-incoming heavy forward
rates (without final heavy Pauli) also stay within .2% of MB
on the refined grid at these three states. The heavier
scalar-incoming fold has more residual quadrature dependence
than the exact thermal ratio, which is why no tighter heavy
accuracy is claimed. This is a pointwise rate audit, not a
quantum integration of a time-dependent distribution or proof
that these factors remain small in a more populated reservoir.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_thermal_occupation_rate_cross_reference.py \
 results/warm-prompt-reheating-expanding-quantum-sensitivity-2026-11-26.json
~~~

[Independent brackets and quadrature controls](results/warm-prompt-reheating-thermal-occupation-rate-cross-2026-11-26.json).
