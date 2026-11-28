# Tagged momentum dependence of finite parent relaxation

Repeat the [frozen-bath real-parent test](INDEPENDENT-EARLY-POLE-SCALAR-PARENT-TRANSIENT-2026-11-27.md) at 18 tagged momenta between 1.1e−6 and 9.5e−6, with the largest previously studied exposure, 7.2e−5. The archived charge and scalar populations are fixed. Only the explicit parent can change by the one reciprocal decay and inverse-decay term. The matched control keeps even the parent frozen.

| Archive | Lowest initial / final net rate over the sampled momenta | Range of accumulated change relative to a frozen parent |
| --- | ---: | ---: |
| Fine | +3.73410 / +3.71437 | −0.563% to −0.265% |
| Third | +3.56743 / +3.52644 | −0.576% to −0.201% |

The largest reduction is at the lowest tagged momentum in each archive. The net remains positive initially and after the modeled relaxation at all sampled momenta. Refining the quadrature order 24 to 32 changes the reduction by at most .00076 percentage points. This is a pointwise momentum scan, with no integrated abundance, no assertion between samples, and no dynamically induced changes in the other populations. The nearly stationary-parent source from the pole identity is not applicable over this finite interval.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
 for n in 24 32; do
  python -W error early_pole_archive_scalar_pole_factorization_reference.py "$g" --order "$n" \
    --q 1.1e-6 1.5e-6 2e-6 2.5e-6 3e-6 3.5e-6 4e-6 4.5e-6 5e-6 5.5e-6 6e-6 \
        6.5e-6 7e-6 7.5e-6 8e-6 8.5e-6 9e-6 9.5e-6 --exposure 0 7.2e-5
 done
done
```
Raw: `results/early-pole-archive-scalar-parent-momentum-scan-2026-11-27.json.gz`.
