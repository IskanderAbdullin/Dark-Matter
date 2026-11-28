# Same dark counts and energy do not fix the heavy regeneration rate

On the [prompt source spectra](INDEPENDENT-WARM-REHEATING-HEAVY-GAIN-PROBE-2026-11-26.md),
test an intentionally different spectral closure. At
each labeled visible temperature, construct the
maximum-entropy Bose and Fermi populations in a
common momentum box at exactly the same separate
scalar and stable numbers and the same summed
dark energy. There is one common dark temperature
and independent chemical potentials. These are
parameters of a comparison, not a measured cooling
history or evidence that equilibration has occurred.

Evaluate the same three Pauli-free regeneration
channels, rather than fixing the cross section
at a typical momentum. Repeat using two momentum
grids and angle rules. Matching the total energy
removes the trivial change from adding heat; it does
not fix the species energy split or the population
above the invariant production threshold.

| Reheat → visible T (MeV) | Reference dark T / visible T | Gross count gain divided by the original passive gain |
|---:|---:|---:|
| 50 → 35 | 1.0102 | 0.9300 |
| 100 → 85 | 0.9150 | 1.3064 |
| 100 → 75 | 0.9112 | 1.3120 |
| 100 → 60 | 0.9042 | 1.2616 |
| 100 → 50 | 0.8961 | 1.1594 |
| 100 → 35 | 0.8704 | 0.8314 |

The reference shape does not quench the fast
gross rate at the hotter epochs despite its
common temperature being less than the visible
temperature. At 60 MeV it gives .006461 units of
extra Q per entropy per Hubble time instead of
.005122. Its local signed scalar energy transfer
per Hubble time is 19.5 times its scalar energy;
the original passive value is 45.4 in its own
energy normalization. Both are instantaneous
rate diagnostics and both require evolution of
the spectrum in any full calculation. This
does not bound the signed net rate with all inverse
and competing channels, or say which of these
two distributions is realized.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error warm_prompt_reheating_heavy_gain_probe_reference.py EOS.late.dat local-spectra > passive.json
python -W error warm_prompt_reheating_matched_entropy_gain_reference.py EOS.late.dat local-spectra passive.json
~~~

[Energy closure, composition and rate controls](results/warm-prompt-reheating-matched-entropy-gain-2026-11-26.json).
