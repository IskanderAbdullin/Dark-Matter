# Positive elastic entropy and reduced total entropy on a different trajectory

Reintegrate the direct combined history to log(a/a_start)=0.25 and retain separate comoving timelike and elastic entropy integrals. Keep the initial state, EOS, number and heat kernels identical to the accepted direct benchmark. Independent BDF and Radau/phase-space resolutions agree on the channel integrals within 8.60e-15 (normalized to initial entropy); summing the channels reproduces the accepted total.

At the terminal clock, the extra elastic term contributes +2.34429e-9 of the initial entropy. The timelike term accumulated on the combined temperature/chemical trajectory is 2.64924e-9 smaller than on the separate timelike-only trajectory. Thus the elastic contribution offsets about 88.5% of the reduction, leaving a total difference of -3.04949e-10. Both individual integrals are positive; only the difference between two distinct histories is negative.

```sh
cd kinetics
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_direct_channel_entropy_expansion.py \
 eos-range.dat results/fd-portal-direct-combined-expansion-tracking-2026-11-28.json
```

The EOS is the same locally normalized authors' table as in the closed-bath benchmark. This result describes the controlled vacuum charged-lepton history, not the full entropy budget of the universe.
