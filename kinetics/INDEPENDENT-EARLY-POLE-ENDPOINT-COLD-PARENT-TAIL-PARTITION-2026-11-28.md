# The extreme retained parent tail is important on individual rays, not necessarily after angular integration

For the independently evolved fine and third four-charge-subcell states at x=.05, the [cold kinematic threshold](INDEPENDENT-EARLY-POLE-ENDPOINT-PAIRED-FOUR-CHARGE-LADDERS-2026-11-28.md) depends strongly on the scalar–photon angle. Keep each state's optical frequency, residue, scalar mass and reconstructed charged occupations fixed at k=.0015. Partition the existing radial integral by parent energy without altering or extending occupation support. Add each requested energy boundary to both the radial panels and the exact polar opening knots. The same spinor evaluation then accumulates the three affine coefficients by parent-energy band. Compare the sum of those bands against the independent unsplit kernel and refine both polar and radial orders; the largest reported relative drift is 1.01e−12 (scaled to the corresponding total coefficient), with independent unsplit closure ≤1.71e−13.

The numbers below are fractions of the **retained forward coefficient** contributed by E_parent/T≥28. An angle-averaged coefficient integrates all allowed cosines at the given q. A single ray estimates neither that angular integral nor the physical forward rate.

| Card | q | Full cosine average | Fixed cosine .9 |
| --- | ---: | ---: | ---: |
| Fine d4 | 3.5e−6 | .0731% | .3071% |
| Third d4 | 3.5e−6 | .2523% | 1.2801% |
| Fine d4 | 1.5e−5 | .0060% | 30.0593% |
| Third d4 | 1.5e−5 | .0144% | 100% |

On the last two fixed rays at q=1.5e−5, the archived support admits a forward coefficient with exactly zero inverse coefficient at the verified orders. The zero reflects the specified retained occupations and is not a physical absence of the inverse reaction. The receipt uses `null` for its undefined inverse fraction. The angle-averaged inverse is nonzero on both cards.

The calculation shows why a high cutoff sensitivity on individual rays does not by itself measure sensitivity of the complete cosine average in these frozen states. “100%” means all forward contribution *in the declared archive*, not that the physical process is dominated by this band. The computation provides no bound beyond the original support ceiling, across all q or k, or for a self-consistent optical solution.

From `kinetics/`, reproduce the [receipt](results/early-pole-endpoint-cold-parent-tail-partition-2026-11-28.json) with both saved states:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  early_pole_endpoint_cold_parent_tail_partition.py \
  results/early-pole-cell-endpoint-{fine,third}-d4-x05-state-2026-11-28.json.gz
```
The optional energy partition requires coverage of the entire retained support. The default unpartitioned kernel retains the original raw coefficients byte-for-byte in the archived mirror fixture.
