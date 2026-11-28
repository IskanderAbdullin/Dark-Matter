# The local stationary lag progressively overstates the live chemical offset

Join the exact two direct FD histories at log(a/a_start)=0.25, retaining their actual particle-count, energy and entropy trajectory benchmarks. At each independently verified archived temperature, solve the instantaneous algebraic equilibrium-force system without replacing the live histories. This separates a local coefficient check from the actual lag built up over time.

In the combined history, the magnitude of live eta is about 86.6% of the algebraic value at log(a/a_start)=0.05, 53.8% at 0.25, and 16.1% at 0.55. At the last anchor, the live value is -0.03919, while the stationary algebraic calculation gives -0.2434. The latter is outside a small-force approximation. A similar progression occurs in the timelike-only path. The two independently evaluated phase-space orders resolve the instantaneous coefficients at the measured states.

```sh
cd kinetics
d=2026-11-28
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error fd_portal_live_histories_vs_frozen_local_lag.py \
 eos-range.dat "results/fd-portal-direct-combined-expansion-tracking-$d.json" \
 "results/fd-portal-direct-expansion-tracking-continuation-$d.json"
```

Use the prior local EOS normalization. These are finite vacuum charged-lepton calculations, not a frozen-displacement model or a decoupling or abundance inference.
