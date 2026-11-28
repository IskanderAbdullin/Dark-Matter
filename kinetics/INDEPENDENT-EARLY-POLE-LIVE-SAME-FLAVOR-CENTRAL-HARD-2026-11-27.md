# A live central-angle same-flavor control with the other hard channels

Add the finite reversible same-flavor Møller/Bhabha sector to the shared cold decay/inverse, photon pairs, Compton and flavor conversion. The reference does **not** include either soft elastic sector in either branch. The compared difference is the same-flavor vacuum hard proxy only, with s≥1 and both dimensionless t/u transfers ≥.15 in the bare event geometry. This sector therefore does not double-count the earlier isolated small-Q t-pole window: the corresponding bare spatial transfer is at least the square root of .15. It does not fill the intervening angular or momentum support.

| Archive and exposure | Baseline cold gain | With separate same-flavor hard sector | Relative change |
| --- | ---: | ---: | ---: |
| Fine, 2e−6 | 1.57100135e−7 | 1.57190188e−7 | +0.05732% |
| Third, 2e−6 | 1.44931820e−7 | 1.44940003e−7 | +0.00565% |
| Fine, 2e−5 | 1.55558093e−6 | 1.56435675e−6 | +0.56415% |
| Third, 2e−5 | 1.43006958e−6 | 1.43119337e−6 | +0.07858% |

On fine at 2e−6, changing the angular quadrature from (1,3,4,1) to (1,4,5,2) changes the fractional effect to +0.06311%. The angular check supports the positive small correction on this fixed archived state; it is not a full radial or hard-mesh convergence test. The gains do not imply a negligible same-flavor sector in general: an earlier [separately held-screened calculation](INDEPENDENT-EARLY-POLE-WITNESS-ELASTIC-SECTOR-ATTRIBUTION-2026-11-27.md) used a deliberately different input and investigated much smaller transfers. Here the cold gain is cumulative, all hard and scalar moments evolve, and the combined hard count and independent energy/entropy ledgers close. The finite cut, no soft scattering, no expansion and fixed channel choices limit the inference.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
for g in fine third; do
  for T in 2e-6 2e-5; do
    python early_pole_live_same_flavor_hard_joint_reference.py "$g" --duration "$T"
  done
done
python early_pole_live_same_flavor_hard_joint_reference.py fine --level 2
```

Raw values: `results/early-pole-live-same-flavor-hard-joint-2026-11-27.json`.
