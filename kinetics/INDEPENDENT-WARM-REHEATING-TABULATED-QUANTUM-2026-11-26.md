# Re-evaluate the quantum collision factors along the separate 60 MeV seed path

Start from the same imposed 100→60 MeV seed of the
[short expanding comparison](INDEPENDENT-WARM-REHEATING-EXPANDING-SLOW-2026-11-26.md).
This is distinct from the [near-empty early history](INDEPENDENT-WARM-REHEATING-EARLY-SEED-SENSITIVITY-2026-11-26.md).
Evaluate rates on selected checkpoints of the current
thermal path, interpolate their positive forward
multipliers, evolve again and recompute the table.
Use the exact chemical-affinity ratio for each
reverse rate. The light-to-scalar branch directly
folds both FD and BE factors. The heavy forward
branches include the incoming distributions and
a midpoint of the final heavy Pauli multiplier
bounds, not a full heavy blocking integral.

| Path table and quadratures | Visible T at phase factor .28 (MeV) | S, scalar-equivalent count per visible entropy, at that point |
|---|---:|---:|
| 9 calibration states, coarser | 51.7624 | 1.39355e−4 |
| 13 calibration states, finer | 51.7715 | 1.39342e−4 |
| 17 calibration states, finer | 51.7715 | 1.39342e−4 |

The terminal dark temperature in the finest trial
is 15.9404 MeV. The refined terminal forward
multipliers over MB are 0.998913, 1.000712 and
1.009252 for light-to-heavy, scalar-to-heavy and
light-to-scalar respectively. The tabulated factors
are evaluations on a nearby reconstructed path,
not a full quantum Boltzmann solution. These results
do not join the independent early run to this seed.

Three additional held-state diagnostics:
- The fast prompt heat gap stays about 1.81e9 to
  3.45e10 times the visible Hubble rate from the
  first to last stored samples. The finite-difference
  quantum energy response checks the corresponding
  linear conductance. The result concerns only the
  relative heat mode and is not a spectrum of the
  full collision operator.
- Include both dark FD blocking and inverse reactions
  in the leptonic visible inflow. The net energy
  inflow stays within 0.03% of the empty-dark result
  on the sampled late states. A manufactured neutral
  common-temperature state has zero net source.
- Holding total energy and the fast light invariant
  fixed, solve the coupled thermal count with scalar
  chemical potential at its mass, using an explicit
  infrared quadrature. At the imposed .28 stops the
  actual invariant count is 48.85% (early path)
  or 49.00% (60 MeV path) of the corresponding
  saturated thermal normal-phase capacity. This
  excludes interpreting the .28 cutoff itself as
  saturation. It does not exclude nonthermal
  condensation or show how either path evolves later.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error warm_prompt_reheating_tabulated_quantum_reference.py EOS.late.dat local-spectra/0.100-0.060-171.json > path.json
python -W error warm_prompt_reheating_quantum_heat_gap_reference.py EOS.late.dat path.json
python -W error warm_prompt_reheating_visible_blocking_reference.py path.json
python -W error warm_prompt_reheating_saturation_distance_reference.py EOS.late.dat path.json
~~~

[Tabulated path and refinement](results/warm-prompt-reheating-tabulated-quantum-2026-11-26.json);
[Bose–Fermi heat mode](results/warm-prompt-reheating-tabulated-heat-2026-11-26.json);
[held source and saturated-capacity checks](results/warm-prompt-reheating-held-source-and-saturation-2026-11-26.json).

The saturated scalar quadrature also reproduces
the leading small-mass corrections to the independent
massless Bose count and energy integrals on a separate
check. These small-mass expansions are not substituted
for the finite-mass saturated-state solution.
