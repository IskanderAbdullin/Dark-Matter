# Tagged early scalar Yukawa elasticity below the parent pole

Fold the positive tree-level scalar–daughter Yukawa elastic channel on the two early archived nonthermal distributions. Both fermion legs are confined to p∈[.45,5]. The invariant guard sqrt(s)/M is .70, .82 or .90, explicitly excluding the real parent pole already represented by decay and inverse decay. This is a frozen, tagged scalar calculation, not an evolved transfer between populations. Momenta, rates and temperatures follow the archived dimensionless conventions. Both scalar partner domains (below or above 1e−5) are separately retained.

The finite-q reduction uses exact center-of-mass boost and azimuthal geometry, splitting the energy panels at both outgoing occupation step boundaries. At q=0 it recovers the independent rest reduction for emission, loss and their normalized difference. The invariant spin sum matches an independent explicit Dirac trace to 3e−15; the common-temperature FD/BE packets close pointwise. Gain minus loss is evaluated as a Pauli difference and two remaining seed terms, so almost-equal large stimulated weights are not subtracted directly. The scalar rest occupation is extrapolated from the first archived bin only for the exact zero-momentum control.

At the widest .90 guard, the matched vacuum-width normalization gives:

| Archive | Tagged q | Net from partners below 1e−5 / Href | Net from partners above 1e−5 / Href | Direct decay/inverse per tagged occupation / Href |
| --- | ---: | ---: | ---: | ---: |
| Fine | 0 | +9.08e−18 | −3.02e−18 | +4.64 |
| Third | 0 | +9.66e−18 | −3.04e−18 | +3.57 |
| Fine | 4.1e−6 | −1.69e−25 | −3.02e−18 | +4.42 |
| Third | 4.1e−6 | −1.80e−26 | −3.04e−18 | +3.79 |

The additional q=2e−6 and 8e−6 controls and all smaller guards keep the net off-pole term negligible compared with the direct source. The largest individual normalized elastic gain or loss in these panels is below 4.1e−10 per reference H. Changing quadrature from 5 to 8 changes the dominant exterior rate by less than 3e−5 relative; the tiny interior remainder is much less well converged relative to its own magnitude and is not treated as a measured transport coefficient. The net test is restricted to the stated vacuum-width normalization, held occupations and cuts. It does not constrain unresolved near-pole elastic transport, virtual and in-medium effects, or a coupled kinetic history.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
  for n in 5 8; do
    python early_pole_archive_scalar_below_pole_finite_q_reference.py "$g" \
      --order "$n" --q 0 2e-6 4.1e-6 8e-6
  done
done
```

Raw: `results/early-pole-archive-scalar-below-pole-2026-11-27.json.gz`.
