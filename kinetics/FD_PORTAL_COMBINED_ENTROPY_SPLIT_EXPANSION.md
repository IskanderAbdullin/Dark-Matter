# Combined timelike and elastic FD collisions in a short expanding split

Extend the [timelike split](FD_PORTAL_TABULATED_ENTROPY_SPLIT_EXPANSION.md) by evaluating direct charged-lepton elastic heat together with the timelike number and heat rates in the energy-conserving collision substep. The nonsplit differential reference independently evaluates both channels at the two actual bath temperatures. Both models use the same tabulated visible EOS and the same pressure-work and entropy budgets. The generalized code reproduces the earlier timelike-only output byte for byte.

Across the three expanding checks, time-step error halving was consistent with second order (ratios 4.0000–4.0735). Both integrated budgets close with errors below 1e-11 on the sampled steps; the scale-factor interval is at most log(a/a_start)=0.01187. With the same physical duration at 0.30 GeV, elastic heat increases the produced comoving entropy divided by initial entropy by 1.295e-6, and changes the terminal dark log-temperature by -6.822e-5. No same-duration comparison is inferred from the 0.68825-GeV pair because the collision-rate time caps differ.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest -q test_fd_portal_tabulated_entropy_split_expansion
python -W error fd_portal_combined_entropy_split_expansion.py eos-range.dat \
  --timelike results/fd-portal-tabulated-entropy-split-expansion-2026-11-28.json
```

The EOS input is documented in the [closed tabulated-bath benchmark](FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md). This does not include all physical collision channels, model the full expansion, or determine a relic density.
