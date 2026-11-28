# The verified fine past-support bridge changes the conditional roots by orders of magnitude

Keep the original archived fine endpoint, scalar mass 0.00087205594, and transverse pole at k=.0015 unchanged. Directly integrate the previously missing adjacent band q=[.0003,.0015]. It contains 108 past-directed panels on both angular/radial meshes. All affine and centered scalar-energy ledgers close; the maximum relative coefficient drift across meshes is 8.06e−6 and the independent charged-energy ledger drifts by at most 6.52e−6. The original cold block is exactly the one in the previously accepted low-range receipt.

The contiguous union q=[1e−6,.0015] is now available for the fine chart. Conditional neutrality means the fixed photon occupation that makes the chosen instantaneous source vanish; it is not the physical occupation. Neither earlier frozen proxy was an all-momentum prediction.

| Held source made neutral | Old q=[1e−6,.0003] | Joined q=[1e−6,.0015] |
| --- | ---: | ---: |
| Scalar number | 1.57912e−8 | 0.20452968 |
| Scalar energy | 1.62340e−8 | 0.20589206 |
| Charged energy | 1.55171e−8 | 0.20204606 |

At the joined number-neutral root, the photon-energy source is zero to roundoff, but scalar and charged energies still exchange equal and opposite amounts, about 7.71e−10 in this fixture. The adjacent bridge accounts for more than 99.86% of integrated direct number response on this union, and nearly all its inverse response. This explains why small-band neutrality is not robust to adding the measured support. It does not infer any physical abundance, evolve the photon bath, determine the third-grid answer, or fill the q>.0015 gap. Independently available higher intervals are deliberately not joined through missing support here.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results
cold="$p/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json"
low="$p/finite-mass-qed-fixed-pole-cold-to-shoulder-gap-$d.json"
python -W error finite_mass_qed_native_scalar_bridge_receipt.py past-bridge \
 --root "$p" --cold "$cold" --meshes 0 1 --grids fine >/tmp/fine-past-bridge.json
python -W error finite_mass_qed_native_past_bridge_partial_join.py \
 --low "$low" --cold "$cold" --bridge /tmp/fine-past-bridge.json
```
