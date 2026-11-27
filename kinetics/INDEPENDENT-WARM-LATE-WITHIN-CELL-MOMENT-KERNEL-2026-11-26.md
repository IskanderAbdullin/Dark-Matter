# Within-cell shape changes the current with the same finite-volume moments

At the archived .0025 GeV half-width frozen root of the
[hard-packet scan](INDEPENDENT-WARM-LATE-TAIL-SPECIES-EXCHANGE-2026-11-26.md),
change a single species inside only one original occupied
cell. Leave the other species and every other original cell
unchanged. Place three adjacent tiles inside
the chosen cell strictly above the physical cut C.
On each tile multiply the initial occupation by a separate
constant; the multiplier is one away from the tiles.

Choose the three coefficients from the null vector of the
number and (E−C) moments of those tiles. Check both with higher
moment quadrature: each residual is less than 5e−16 of the
absolute weighted cancellations. As the affected support is
all above C, the perturbed states preserve each original
cell's full number and energy, and also preserve the
above-cut count, full energy, and convex excess of each
species. The perturbed scalar occupations remain positive
and the perturbed fermion occupations remain between zero
and one. These positive piecewise occupations are valid
frozen functions but need not be the entropy-maximizing
shape for those moments.

| Changed species | Amplitude −0.5 signed current | Amplitude +0.5 signed current |
| --- | ---: | ---: |
| Scalar | +1.92063e−13 | −1.92063e−13 |
| Light fermion | +2.59176e−13 | −2.59176e−13 |

The four signed rates are in units of initial dark energy
times H, computed by splitting the exact angular integrals
at the new tile walls. Three angular quadratures agree and
the zero-amplitude tiling recovers the unaltered result.
The unnormalized comparison is close to a cancellation:
the baseline signed rate is near −1e−19, while the sum of
gross growth and debit rates is about 5.47e−7 in the same
units. This shows that **the original two moments in every
cell still do not specify the sign near a frozen cancellation**,
not that the choice of within-cell shape creates a large
net energy transfer. Refining the representation or adding
closure information can distinguish the perturbations.
This does not test whether these shapes are produced by the
physical time evolution, and makes no structure-formation
claim.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
python -W error warm_late_within_cell_excess_current_reference.py EOS.dat \
  results/warm-corrected-anchor-signed-floor-2026-11-26.json \
  results/warm-multishell-corrected-q64-fine-stage-48-2026-11-26.json.gz \
  results/warm-late-frozen-hard-center-roots-2026-11-26.json
~~~

[Tile coefficients, null-moment audit, and signed/gross rates](results/warm-late-within-cell-excess-current-2026-11-26.json).
