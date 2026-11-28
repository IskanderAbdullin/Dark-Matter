# Independently re-integrated factors distinguish aperture and scalar population at four cold points

Match the two approved native D2 endpoints and their complete raw panel receipts before making crossed evaluations. At each fixed q re-integrate all eight past-directed gross bare states, varying the scalar gap, photon frequency/residue, and charged occupation closure. Apply the explicit scalar factors from the separately stored scalar populations. Average every log replacement order. Independent direct calculations of the fully weighted native endpoints recover all four forward and inverse ratios.

| q | Third / fine forward response | Third / fine inverse response | Scalar-gap share of bare forward log change | Explicit scalar share of total forward log change |
| ---: | ---: | ---: | ---: | ---: |
| 2.5e−6 | 0.19512 | 0.16347 | 80.32% | −8.75% |
| 3.5e−6 | 0.29909 | 0.25127 | 80.34% | −44.93% |
| 4.7e−6 | 0.03298 | 0.02783 | 80.32% | +49.87% |
| 8.2e−6 | 0.75373 | 0.64303 | 80.28% | −459.98% |

A negative share means that an increase in the explicit factor offsets part of a decrease in bare response. Values with magnitude over 100% can arise from cancellation in the total log change. In this test the mass-gap contribution dominates the **bare** response at all four points, consistent with the independent collinear opening check. The explicitly retained scalar profile is not a uniform correction: at 4.7e−6 it increases the suppression, at the other points it offsets part of it. The direct fully weighted ratios agree with the factored corners to 1.47e−13 relative. This comparison is only at four fixed q values; it does not decompose the integrated −83.1% response or imply a physical cosmological result.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results; d=2026-11-28
python -W error early_pole_endpoint_cold_fine_third_point_factors.py \
 --fine "$p/fine-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --third "$p/third-d2-endpoint-continuum-cold"-{1,2}.jsonl.gz \
 --states "$p/early-pole-cell-endpoint"-{fine,third}-d2-x05-state-$d.json.gz \
 --archive "$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json" \
 --cdf "$p/early-pole-cold-fine-third-d2-native-cdf-$d.json"
```
