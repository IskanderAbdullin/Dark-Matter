# Repeat the two-point aperture attribution on more refined evolved cards

Apply the [same frozen sequential replacement](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-CHARGE-REFINEMENT-2026-11-28.md) to the new fine four-charge-subcell and third two-charge-subcell x=.05 states. Starting from each grid's original archive, replace the scalar mass, then the transverse photon pole and residue, then both continuum charged occupations. Rates are evaluated at k=.0015; the two q points are not a band integral. The angular/radial mesh comparison was below the same 2e−6 relative gate on all stages.

| Evolved state | q | Mass-only factor | Subsequent optical factor | Subsequent charge-occupation factor | Total forward factor |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fine, four charge subcells | 3.5e−6 | .07683 | .98960 | 1.0000000021 | .07603 |
| Fine, four charge subcells | 8.2e−6 | .09945 | .99073 | 1.0000000024 | .09853 |
| Third, two charge subcells | 3.5e−6 | 1.67073 | .99687 | 1.0000000019 | 1.66550 |
| Third, two charge subcells | 8.2e−6 | 1.58464 | .99719 | 1.0000000024 | 1.58018 |

In the original fine archive the cosine-zero parent threshold at q=3.5e−6 is E/T=18.35. In the fine four-subcell card it rises to 20.60, consistent with the greatly reduced forward coefficient. The same probe on the one- and two-subcell fine cards gave factors 6.459 and 1.304 relative to that **same** original fine archive. Increasing charge resolution can therefore reverse the earlier inferred direction of this pointwise comparison. It is not justified to use a small time-step defect, an optical envelope at fixed evolved state, or a relatively small change in total yield as an uncertainty for the cold response. Attribution refers to this path only and is not a split of coupled physical plasma uncertainties.

From `kinetics/`, reproduce the [receipt](results/early-pole-endpoint-fine4-third2-cold-aperture-2026-11-28.json) using both archived x=.05 states decompressed to JSON:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  early_pole_endpoint_cold_aperture_attribution.py FINE_D4_X05.json THIRD_D2_X05.json
```
