# Formal photon-only test across the exact emission endpoint

The [cold boosted reference](INDEPENDENT-FINITE-MASS-QED-BOOSTED-HIGHER-BATH-2026-11-27.md)
is an asymptotic expansion. Here integrate the same **formal bare
photon-only** O(g²) comparison without expanding the real phase
space, to see where extrapolating that expansion fails. Photon
absorption has no upper energy cutoff in this two-body channel.
Stimulated emission closes at k/M=[1−(m+s)²/M²]/2. Subtract the
common signed leading soft term before integration; keep the
independently matched external mass and spin response.

Split at the emission endpoint and retain the absorption and
virtual subtraction above it. Check the exact angular primitive
against a separate logarithmic-angle quadrature for signed photon
energies from −240 M (absorption in the signed convention) to
near the emission endpoint. The tested angular relative differences
are below 2.2e−14. Independent adaptive and finite-panel Bose
integrations of the full remainder agree within 7e−9 in the
stated dimensionless rate coefficient.

Quoted coefficients below divide the relative correction by
[g²/(4 pi²)] (pi²/6) (T/M)². For the fine mass fixture in the
parent proper clock with the parent at rest in the bath, the
values at T/M=.04, .1 and .8 are about −2.03024, −2.27744 and
−4.45697. The corresponding coefficients for beta=.8 are
−2.14037, −2.67526 and −4.87000. Where the exact emission
endpoint matters, simply adding the next cold recoil power
does not represent the comparative rate. With m/M=.2,
s/M=.65, a formal high-temperature chiral comparator crosses
zero near T/M=17.38, while the scalar comparator at the same
masses has no sign change on the sampled interval. This is a
counterexample within the simplified formula, not a physical
hot-plasma sign claim.

No equilibrium charged bath, Pauli blocking, thermal scalar
population, resummed photon screening, broadening, or
self-consistent kinetic background is included. At large T
this exercise checks the limits of a bare comparison, not the
sign or magnitude of a hot archive production rate.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_full_support_photon_bath_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-full-support-photon-bath-2026-11-27.json).
