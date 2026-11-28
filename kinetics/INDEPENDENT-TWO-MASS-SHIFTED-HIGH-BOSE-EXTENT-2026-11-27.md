# A nonlinear near-balanced event with sub-ulp fermion changes

To complement the [high-occupation arithmetic](INDEPENDENT-TWO-MASS-SATURATED-BOSE-GRADIENT-2026-11-27.md),
construct a single manufactured Bose/Pauli event with Bose count
(1/epsilon) and two finite Fermi occupations near a common
reference f. Capacities are unity in this toy; the large count
is not inferred from a physical spectrum.

Represent the fermion counts as f+εP and f+εQ rather than
discarding the small corrections in absolute occupations.
For forward reaction extent ξ=εx, the shifted coordinates
become P−x and Q+x, with Bose count (1/epsilon+epsilon x).
Keep the shifted logits in log1p form and evaluate the signed
quantum flux with expm1. Use time s=t/ε, so dx/ds=S.
In the ε→0 control its leading form is

    dx/ds = f(1−f) + P−Q − 2x.

This gives an independent analytic limit for both forward
and inverse starts. The scaled entropy gain is
Δ𝒮/ε²=∫(A/ε)dx. Compare that integral to a separate
120-digit direct subtraction of the two endpoint entropies.

For ε=10⁻²⁴, absolute float64 Fermi occupations cannot resolve
the count changes. The shifted calculation still returns
x=.26611676 and −.13785555 for the two starts, agrees with
the analytic limiting solution to below 5.8e−15 in both, and
closes the scaled entropy comparison within 3e−16 relative
across the tested ε range. No information lost in the absolute
occupation representation is reconstructed: the shifted
perturbation has to be present in the input.

This is a one-event stress test of numerical coordinates,
not an evolved high-occupation spectrum or a full nonlinear
solver with many competing channels.

Replay from the kinetics directory with NumPy, SciPy and mpmath:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_shifted_high_bose_extent_reference.py
```
[Full receipt](results/two-mass-shifted-high-bose-extent-2026-11-27.json).
