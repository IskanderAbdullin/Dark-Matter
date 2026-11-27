# Lagged updating of the isolated scalar scattering proxy preserves the tested signs

Extend the [live reciprocal decay and redistribution test](INDEPENDENT-EARLY-POLE-JOINED-SCALAR-SELF-DECAY-PROXY-2026-11-27.md). Before each scalar redistribution half step, rebuild the screened scalar vertex and its scattering dispersion using the current scalar occupations and Hartree mass; hold them within that half step. The hard states remain live through the reciprocal decay exactly as before. This checks sensitivity to refreshing the proxy, not a controlled physical collision amplitude. The no-redistribution trajectories are exactly identical to the preceding coarse runs.

| Witness at 2e−7, Q16 and 8 macro steps | No scalar redistribution | Held proxy | Lagged refreshed proxy |
| --- | ---: | ---: | ---: |
| Fine 97% parent | −1.36436e−9 | −9.88457e−10 | −1.01183e−9 |
| Third 97% parent | +1.30336e−8 | +1.04883e−8 | +1.07234e−8 |

The refreshed proxy still reduces the magnitude of the fine witness loss by 25.8% and the third witness gain by 17.7% relative to no scalar redistribution. Refreshing changes the accumulated witness responses by about 2.4% and 2.2% relative to holding the initial proxy. The final scalar shape L1, normalized to the initial count, changes from .859 to .792 in fine and from .705 to .633 in third. Doubling the steps and increasing angular order together gives −1.01675e−9 and +1.07094e−8 in the refreshed branches. Both signs survive. This joint refinement does not separately quantify time and angular errors.

The work against the evolving Hartree gradient and the independent entropy production remain tracked. The exact endpoint gradient serves as a separate audit including the fixed spectators, not as an imposed change to the scattering geometry. The broad reasons for model uncertainty remain: only the finite cold scalar support, the illustrative scalar vertex, and reciprocal decay are included. No screened charged scattering or FRW occurs in these comparisons. A result for the model with the additional charged processes requires a separate matched calculation.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
for g in fine third; do
  python early_pole_joined_scalar_self_proxy_reference.py "$g" --refresh
  python early_pole_joined_scalar_self_proxy_reference.py "$g" --refresh --steps 16 --order 10
done
```

Raw: `results/early-pole-joined-scalar-lagged-refresh-2026-11-27.json.gz`.
