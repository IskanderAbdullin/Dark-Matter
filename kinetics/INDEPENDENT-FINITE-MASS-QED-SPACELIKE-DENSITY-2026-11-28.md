# Free spacelike density statistics and the longitudinal Ward block

This is the charge-density counterpart of the [transverse source](INDEPENDENT-FINITE-MASS-QED-TRANSVERSE-STATISTICS-2026-11-27.md), not the population or residue of a longitudinal collective pole. Hold the same two archived massive isotropic charged species and mirror their antiparticles as before. Write K = (ω, 0, 0, k), u = ω/k, s = k²−ω², a = 1−u², b = k√(1+4m²/s), and E₀ = (b−ω)/2. With both charged legs on shell the density trace is

```text
M00(E) = 8E(E+ω) − 2s
       = 8[m²/a + (E−E₀)(E−E₀+b)].
```

The same energy interval and one-dimensional phase measure as for the transverse cut apply. Current conservation fixes the 03 and 33 entries to u and u² times the 00 entry (contravariant Cartesian currents). Explicit spinors verify these identities as well as the trace and null Ward contraction to 8.89e−16 relative.

Resolve each forward/inverse Fermi boundary in energy separately. The rearranged dissipation follows from M00(E)−M00(E−ω) = 16Eω:

```text
D00/ω = integral[E₀,E₀+ω] M00(E) f(E) dE / ω
        + 16 integral[E₀+ω,∞] E f(E) dE.
```

Hence the free isotropic density cut has nonnegative damping, even for a nonmonotone occupation. Its first term survives at fixed u as k→0 if m≠0: after division by 8 the soft denominator is E*² f(E*) + 2∫ E f(E)dE above E* = m/√a. The histogram jump-sum independently gives the same result; if the threshold lies exactly at an edge, half the two one-sided occupations enter because the recoil layer is centered on the threshold at leading order. The numerator after division by 8 is ∫ E² f(1−f)dE above E*. The synthetic edge test explicitly plants a jump at threshold.

The archive swept-weight denominator agrees with direct signed integration to 7.69e−11; its soft limit and the jump check agree to 3.34e−16 on the actual archives. The independently integrated common-temperature Fermi control reproduces the finite-frequency thermal fluctuation law to 6.78e−15. In the test range k≤.01, finite recoil changes the density statistical ratio by at most 8.14e−5 (fine) and 1.007e−4 (third). The soft density ratios near u=.9996 are .94529 and .94950, larger than the corresponding transverse ratios .92585 and .92787 on the same charged occupations. This difference is not a temperature split between equilibrated photon modes.

Run `python -W error finite_mass_qed_spacelike_density_statistics_reference.py` from `kinetics/`. The reproducible receipt is [under results](results/finite-mass-qed-spacelike-density-statistics-2026-11-28.json). No anisotropic instability, longitudinal collective occupation, resummed vertex, interacting width, or matched scalar source is inferred.
