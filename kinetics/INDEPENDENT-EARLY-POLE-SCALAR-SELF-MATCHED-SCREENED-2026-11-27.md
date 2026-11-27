# The first matched screened fine test retains the prepared-witness sign reversal

Add photon pairs, hard Compton, hard flavor conversion, and both finite screened elastic sectors to the live fine 97% witness. Independently vary only the scalar self-redistribution proxy, keeping its vertex and geometry held at initialization. The tested interval is 2e−8, with eight macro steps, three hard subdivisions and sixteen scalar bins. These are the same numerical controls as the earlier [five-channel witness test](INDEPENDENT-EARLY-POLE-WITNESS-ELASTIC-SECTOR-ATTRIBUTION-2026-11-27.md), with added matched scalar on/off branches.

| Matched setting | Cold gain | First positive rate at 2.5e−9 | Final instantaneous rate |
| --- | ---: | ---: | ---: |
| No scalar self | +9.69149e−10 | +.0302631 | +.06067869 |
| Held scalar self added | +9.69880e−10 | +.0303128 | +.06068835 |

The scalar shape L1 reaches .1523, yet the integrated change between the two branches is only +0.0755% in this short control. The small difference is not a converged estimate of a physical correction. In the separate longer 2e−7 control **without either screened elastic sector**, with shared hard pairs, Compton and conversion, the scalar proxy changes the prepared loss from −1.35860e−9 to −9.84003e−10. These windows differ and must not be treated as a matched estimate of the screened suppression of a single term.

Active unlike screening replaces graph weights. Each comparison arm therefore restores the pristine input graph before constructing its response object. Reversing branch order in the small independent fixture gives bitwise identical, sorted outputs. On the screened comparison, total energy/work closure is below 5e−16 and entropy mismatch below 5e−13. An optional per-bin diagnostic records the arithmetic sum of tagged scalar source terms and the hard counts. This test neither establishes full angular support nor constrains late or cosmological behavior. It is a prepared input, and the scalar vertex remains an illustrative proxy.

From `kinetics/`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
python early_pole_scalar_self_branch_order_reference.py
python early_pole_joined_scalar_self_proxy_reference.py fine \
  --duration 2e-8 --hard-slow --screened --mixes .97
python early_pole_joined_scalar_self_proxy_reference.py fine --hard-slow
```

Raw: `results/early-pole-scalar-self-matched-screened-2026-11-27.json.gz`.
