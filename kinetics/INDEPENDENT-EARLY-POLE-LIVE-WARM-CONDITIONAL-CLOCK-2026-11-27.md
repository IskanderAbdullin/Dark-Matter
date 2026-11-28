# Initial density explains the resolved warm-gain normalization change

Use the [short live warm-marker run](INDEPENDENT-EARLY-POLE-LIVE-WARM-REDSHIFT-2026-11-27.md) from the accepted fine final n216 hard snapshot with the same prescribed log(a)=t, six steps, exposure 2e−6, and held polarization. Compare the original radiation-only normalization with dividing the collision rates by the conditional factor (\sqrt{1+\rho_\mathrm{model}/\rho_\mathrm{vis}}), where (ho_\mathrm{vis}=\pi^2 106.75/30) in the existing comoving units. Use the prepared density once or refresh the chosen finite-state density during evolution. The inverse factor multiplies all joined collision terms, not the explicit redshift trace integral. This distinguishes a clock-normalization test within a fixed model from solving cosmological expansion.

| Clock choice | Warm gain |
| --- | ---: |
| Radiation only | 5.49688961028e−13 |
| Conditional initial density held | 5.36094457988e−13 |
| Conditional live density | 5.36094457983e−13 |

The common initial inverse clock factor is 0.97526974019688. Its fixed normalization accounts for a 2.4731% reduction of the warm gain. The difference between refreshing the density and holding it is only about 5e−24 in the warm gain, below what is useful to interpret as a physical effect; the terminal live factor differs from its initial value by about 1.03e−11. All runs pass the original terminal comoving energy, entropy and total hard count checks. No measurable time-dependent clock correction, calibrated cosmology, or late relic prediction follows from this short test.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_frw_reference.py
for k in radiation initial-density live-density; do
  python -W error "$p" "$a" --expands --clock-mode "$k"
done
```

Raw: `results/early-pole-live-warm-conditional-clock-2026-11-27.json`.
