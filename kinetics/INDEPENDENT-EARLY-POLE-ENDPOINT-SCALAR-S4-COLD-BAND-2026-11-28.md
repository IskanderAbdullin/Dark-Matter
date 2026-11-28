# The matched S4 cold band resolves, but its local shape still changes

At x=.05 use the native fine charge-d2 partition in all three evolved spectra. Independently vary the scalar resolution below q=1e−5 (S1/S2/S4), reconstruct the optical pole at k=.0015, and integrate q=[1e−6,1.5e−5]. Both full S4 meshes have 71 contiguous panels, exact state fingerprint and mirror azimuth. The largest relative difference over nine independently integrated number, energy and excess coefficients is 1.52e−13. All three archived topology checks reproduce each selected state and exact panel edges.

| Scalar chart | Panels | Gross forward vs original fine archive | Conditional number-neutral occupation | Relative energy-minus-number split |
| --- | ---: | ---: | ---: | ---: |
| S1 | 49 | 1.30271123 | 2.90893020e−9 | +5.34301e−8 |
| S2 | 57 | 1.37189081 | 2.91523004e−9 | +3.86604e−8 |
| S4 | 71 | 1.36736113 | 2.91357446e−9 | +3.47222e−8 |

The S1→S2 gross forward integral changes +5.3104%, then S2→S4 changes −0.33018%. S2→S4 inverse changes −0.38678%. Nearly unchanged integrated coefficients do not show that the conditional source is pointwise stable.

A separate held-q factorial swaps gap, photon energy/residue, the fitted charged populations, and the scalar occupation multiplier, considering all three-factor orders and separating the scalar factor exactly. At q=3.5e−6 the weighted gross forward change S2→S4 is **+17.484%**, whereas at q=4.7e−6 it is **−41.490%**. The mass-only factors are .9999173 and .9999184 and the photon and charged shifts are smaller, so on these points the within-band occupation is the dominant change. S1→S2 similarly moves these two points in opposite directions; the S1→S2 mass factors are ~1.078, much smaller than the occupation factors (1.743 versus .404). These two samples do not explain the whole band integral; they show how local changes can remain large even when its weighted total is similar. All crossed point quadratures agree within 4.81e−13 relative.

The positive neutrality splits define conditional roots at the retained pole and band, not a predicted photon occupation. These are independent numerical diagnostics, not a CLASS result or an out-of-band/free-streaming prediction.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
p=results/early-pole-fine-charge2-scalar
f=results/early-pole-cell-endpoint-fine-d2-x05-state-$d.json.gz
python -W error early_pole_endpoint_cold_scalar_band_ladder.py \
 --s1 results/fine-d2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s2 results/fine-charge2-scalar2-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --s4 results/fine-charge2-scalar4-endpoint-continuum-cold-{1,2}.jsonl.gz \
 --states "$f" "$p"{2,4}-x05-state-$d.json.gz \
 --archive results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json
python -W error early_pole_endpoint_scalar_band_point_factorial.py \
 "$f" "$p"{2,4}-x05-state-$d.json.gz
```
Rerun a raw S4 panel stream with `early_pole_endpoint_continuum_cold_band_pilot.py` on the uncompressed S4 state, mesh 1 or 2 and `--mirror-azimuth`.
