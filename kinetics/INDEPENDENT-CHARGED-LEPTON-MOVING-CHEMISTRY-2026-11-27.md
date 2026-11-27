# Crossing a prescribed chemical clock in time at fixed thermodynamics

The [static crossover](INDEPENDENT-CHARGED-LEPTON-FINITE-CHEMICAL-CROSSOVER-2026-11-27.md)
tests eigenvalues, not evolving perturbations. Here use the same
entropy-normalized five-moment generator at fixed common temperature,
fixed enthalpy ratio h_V/h_D=10, and fixed wave number k=.00025 in
elastic heat units. Impose an exponential decrease of the full
timelike process multiplier so that its projected chemical rate
runs from 100 k to .01 k. All elastic and thermodynamic coefficients
stay fixed. The control does not use a cosmological background.

Start in the normalized forward acoustic eigenvector at the beginning.
Integrate the complete varying generator in complex coordinates and
replay it with an independent real-block stiff integrator. The
thermodynamic metric is constant; the squared norm loses exactly
the time-integrated positive collision quadratic form. Separately
evolve the instantaneous passive three-mode generator and reconstruct
its two eliminated fast coordinates at the end.

A further comparison keeps only the exponential of the integrated
instantaneous acoustic eigenvalue. This is a diagnostic that omits
changing eigenvector geometry and mode conversion. Compare both the
full-vector norm and the final coefficient of the forward acoustic
mode in the instantaneous eigenbasis; that coefficient is *not* a
standalone conserved energy fraction.

At .015 GeV, the fast e-fold (kτ=.03) leaves the total norm above
.9999995, while the forward amplitude is .9996453 and the
instantaneous eigenvalue-only estimate .9999811. For kτ=30 the
forward amplitude and the instantaneous estimate are much closer
(.9812930 versus .9812928). The passive three-mode model tracks the
full vector in both controls. Across all eight tests, the endpoint
loss ledger closes to 1.54e-8 and the two integrators agree to 7.7e-9. This separates mode rearrangement
from true dissipative norm loss and does not identify a cosmological
transition, FRW evolution, or transferred dark matter power.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error \
  charged_lepton_moving_chemistry_reference.py \
  results/charged-lepton-portal-combined-momentum-2026-11-27.json
```
[Full receipt](results/charged-lepton-moving-chemistry-2026-11-27.json).
