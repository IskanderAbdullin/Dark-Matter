# Direct FD continuation after the earlier fast-clock labels

The exact two direct histories now reach log(a/a_start)=1.20 at visible temperature about 0.29906 GeV. From the x=.96 checkpoint, the additional dark comoving count is 0.8678% (timelike only) or 0.8655% (with elastic scattering) of the original count. At x=1.20 the instantaneous count changes are still 0.03340 and 0.03330 of the initial comoving count per unit log(a/a_start). These are model histories and normalizations, not abundance constraints. The dark heat column is total collision energy transferred to the dark bath, including energy carried by number-changing events.

The minimum of the dark reduced chemical potential is at log(a/a_start)=0.9601145 (timelike) or 0.9905445 (with elastic scattering), with eta=-0.0759358 or -0.0758983. At these minima the count source remains positive, about 0.04018 or 0.03863, and the local entropy production is positive. In the combined case the derivative contributions from adiabatic expansion, number exchange and energy exchange are +0.0019192, +0.0946661, and -0.0965853; their cancellation explains the minimum.

The independent intensive histories agree in their state variables to 9.55e-11; the separately integrated extensive state agrees to 1.64e-11. The two methods locate the minima to 2.56e-09 in their most discrepant reported coordinate. The earlier exact entropy, pressure-work and number budgets are retained.

Scope: the local fast rates discussed earlier were calculated by linearizing near equal temperatures and zero chemical potential. The independent intensive and extensive equations used here instead follow the actual chemical potential and temperature difference. Nonzero flux beyond a local linearized clock label, or at the minimum of eta, is not a cosmological freeze-out criterion. The check excludes hadronic reaction kernels, additional channels, and inference of the final relic abundance.

Numerical method: integrate in log(a/a_start) from the exact accepted .96 endpoint. Compare two quadrature/solver choices (BDF and Radau) and an extensive formulation based on independently inverted comoving dark number and energy plus total energy. Use the original entropy, pressure-work and number ledger gates. Around the eta minimum, reconstruct a dense solution on the verified sign bracket and decompose the actual eta derivative using the thermodynamic response matrix, with expansion, number exchange and energy exchange as separate columns. A zero sum does not set any individual column to zero.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; e=eos-range.dat; p=results
before="$p/fd-portal-direct-expansion-055-096-$d.json"
late="$p/fd-portal-direct-expansion-096-120-$d.json"
python -W error fd_portal_direct_late_tail_local_sources.py "$e" "$before"
python -W error fd_portal_direct_eta_turning_points.py "$e" "$before" "$late"
```

Use the same public EOS normalization documented in `FD_PORTAL_TABULATED_BATH_ENTROPY_AND_MODES.md`.
