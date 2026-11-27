# Occupation factors, not the hard dispersion change, set the initial cold-source contrast

Use the same prepared native and continuous hard states as the [joint moment-matching control](INDEPENDENT-EARLY-POLE-PREPARED-HARD-JOINT-MOMENT-MATCH-2026-11-27.md). Hold the original eight cold scalar nodes, occupations and scalar pole fixed. Cross the variational hard energies used in the packet geometry with the fitted Fermi occupations used for gain and inverse loss, then evaluate the same reciprocal rate integral. The independently evaluated native and continuous rates agree with the archived cold mode sources. A crossed occupation/energy pair is a held-state diagnostic, not a self-consistent state to evolve.

| Hard energies | Hard occupations | Net cold source |
| --- | --- | ---: |
| Native | Native | .07863619025 |
| Continuous | Native | .07863967741 |
| Native | Continuous | .05689797247 |
| Continuous | Continuous | .05690395094 |

Changing only the occupations changes the rate by −.02173822. Changing only the energies at native occupations gives +3.49e−6, with a +2.49e−6 mixed term. The pure dispersion shift is much smaller than the occupied source contrast, and even has the opposite sign. With both inputs changed, the gross gain falls by .00802255 and the gross inverse rate rises by .01370969, accounting together for the net fall of .02173224. Applying only the daughter occupation swap at native dispersion gives net −.01551119, of which .01331208 comes from increased inverse loss. The separate parent swap has net −.00622702, dominated by reduced gross gain.

At this prepared scalar seed the explicit spontaneous rate evaluated at zero scalar occupation is around 1e−16, many orders below the net occupied cold source. That is unlike the [warm-shell quadrature comparison](INDEPENDENT-EARLY-POLE-WARM-SPONTANEOUS-QUADRATURE-2026-11-27.md), where the spontaneous term accounts for almost all the warm marker bias. The seed values and held hybrids in this calculation are model-specific and imply neither a finite-time rate attribution nor a statement about the physical incoming spectrum.

From `kinetics/`:
```sh
OPENBLAS_NUM_THREADS=1 python -W error early_pole_initial_cold_source_factorization_reference.py \
  results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
```

Raw: `results/early-pole-initial-cold-source-factorization-2026-11-27.json`.
