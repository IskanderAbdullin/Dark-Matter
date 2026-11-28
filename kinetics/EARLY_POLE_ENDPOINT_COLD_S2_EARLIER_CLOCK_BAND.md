# Whole-band S2 contrast changes only the earlier .02–.04 clock

Both original and alternate S2 states use the same cold/charge spatial chart and the same 2048-step final segment from x=.04 to .05. Only the preceding .02–.04 interval changes from 2048 to 1024 steps. The alternate final state passes a complete pair of 57 cold quadrature panels per mesh over [1e−6, 1.5e−5]. The maximum whole-band component discrepancy of these meshes is 9.85e−14. The exact state, earlier restart and point metadata bind the checked chronology.

On the complete band, the alternate-minus-original shift splits into a combined bare kernel part and the explicit scalar factor (the latter is (1+f) forward and (f) inverse). The four crossed response profiles are exact panel sums because the scalar chart is unchanged.

| Gross response | Earlier-coarse / original minus one | Explicit scalar share of integrated log shift |
| --- | ---: | ---: |
| Forward number | −0.458954% | 1.11553% |
| Inverse number | −0.498685% | 1.02575% |
| Forward excess | −0.455790% | 1.09618% |

The original and altered forward profiles differ by just 1.413e−5 TV on the common 17 display bins despite the roughly 0.46% shift in amplitude. The independently reconstructed four-corner response agrees across both meshes to 7.16e−14. This earlier timing axis is again mostly a bare-kernel change, unlike the one exactly evaluated spatial replacement order. This test does not split the combined kernel further, bound global time error or establish a physical reaction prediction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/fine-charge2-scalar2
a=results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_continuum_cold_band_receipt.py \
 "$p-prior-clock-endpoint-continuum-cold"-{1,2}.jsonl.gz --archive "$a"
python -W error early_pole_endpoint_cold_prior_clock_population_rectangle.py one --level 2 \
 --original "$p-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --alternate "$p-prior-clock-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states results/early-pole-fine-charge2-scalar2-x05-state-$d.json.gz \
   results/early-pole-cold-s2-prior-clock-x05-state-$d.json.gz \
 --archive "$a" --history results/early-pole-cold-s2-prior-clock-report-$d.json \
 --point results/early-pole-cold-prior-clock-points-$d.json \
 --display results/early-pole-cold-spatial-population-factor-bounds-$d.json
```
