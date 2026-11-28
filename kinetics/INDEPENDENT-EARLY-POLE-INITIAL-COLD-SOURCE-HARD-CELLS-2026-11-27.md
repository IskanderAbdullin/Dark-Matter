# Two retained hard cells per species carry the initial occupation effect

At the frozen native dressed hard energy map, replace the fitted occupation of one hard cell at a time with its continuous value and recalculate the initial cold source. All 20 retained cells in each charged species differ in their fitted logits, but only two per species contribute a detectable substitution response on the tested cold packet support. Keep the other species, the scalar seed and the geometry unchanged. For one substituted species at a time the decay/inverse kernel is affine in its occupation, so the cell changes must sum to the directly evaluated whole-species swap. The maximum numerical additivity defect across mode terms is below 1.4e−17.

| Retained hard momentum cell | Parent occupation swap: net rate change | Daughter occupation swap: net rate change |
| --- | ---: | ---: |
| [1.413889, 1.422222] | −.000470139 | +.001125080 |
| [1.422222, 1.466667] | −.005756885 | −.016636274 |
| Other retained cells | No change detected | No change detected |

The daughter high-cell change comes mainly from a +.01427765 inverse-rate shift; its lower cell partly cancels that suppression. These signed cell contributions agree at radial orders 8 and 12. The support statement applies to the [initial cold occupation effect](INDEPENDENT-EARLY-POLE-INITIAL-COLD-SOURCE-FACTORIZATION-2026-11-27.md) at one held dispersion. It does not mean that the other fitted spectra coincide, that these cells alone control later rates, or that simultaneous parent and daughter changes are generally additive.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for n in 8 12; do
  python -W error early_pole_initial_cold_source_cell_attribution_reference.py \
    results/early-pole-joined-fine-216-memory-2026-11-27.json.gz --order "$n"
done
```

Raw: `results/early-pole-initial-cold-source-hard-cells-2026-11-27.json`.
