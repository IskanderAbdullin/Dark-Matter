# Positive free transverse recoil without subtracting nearly equal rates

The [exact statistical aperture](INDEPENDENT-FINITE-MASS-QED-TRANSVERSE-STATISTICS-2026-11-27.md) verified positive transverse dissipation on two particular isotropic massive archives. For this **free spacelike transverse cut**, nonnegativity is in fact structural for any nonnegative isotropic occupation with convergent moments. The occupation need not decrease with energy.

Let a = 1 − (ω/k)², s = k² − ω² > 0, b = k √(1 + 4m²/s), and E₀ = (b − ω)/2. The spin-summed weight averaged over transverse axes is

```text
M(E) = 2s + 4a (E − E₀) (E − E₀ + b).
```

For E ≥ E₀ + ω the *exact* difference M(E) − M(E − ω) is 8a E ω. Shifting integration variables in the inverse term therefore proves

```text
D/ω = integral[E₀,E₀+ω] M(E) f(E) dE / ω
      + 8a integral[E₀+ω,∞] E f(E) dE.
```

Both terms are nonnegative, and any occupied interval above the aperture opens strictly positive damping. The first term is a thin threshold layer, integrated in coordinates x = (E − E₀)/ω. The second integrates directly over the archived histogram energies. This does not require resolving differences in Fermi factors or taking gain minus loss. Since M(E₀) = 2s vanishes in the soft limit at fixed phase speed, D/ω tends to 8a times the integral of E f(E) dE above the soft energy threshold, or exactly four times the earlier soft D convention.

## Numerical audit

- Both archive endpoints were tested at six phase speeds from .0047 to .9996 and six k values from .0002 to .5. Relative agreement with independently split forward/inverse energy bins is better than 7.7e−11. The larger k values are mathematical checks, not an assertion about the scalar emission aperture.
- At k = 1e−12 the normalized dissipation agrees with the soft limit to 2.52e−13. At k = 1e−22 the symmetric source divided by the stable D/ω still recovers the independent soft statistical ratio.
- Of 288 deliberately nonmonotonic single-cell controls spanning three masses, 236 have positive dissipation and the rest have exact zero support. The direct finite-recoil routine agrees to 6.59e−12 for the sampled non-negligible cases. Random nonmonotonic masks also pass the positive ultrasoft gate.

Reproduce from `kinetics/` with
`python finite_mass_qed_spacelike_positive_dissipation_reference.py`.
The receipt is [under results](results/finite-mass-qed-spacelike-positive-dissipation-2026-11-27.json).

This identity excludes anisotropic states, time-like pair cuts, interacting widths and vertex corrections. It supplies a stable free retarded/statistical normalization, not a pole population or contribution to the matched scalar collision term.
