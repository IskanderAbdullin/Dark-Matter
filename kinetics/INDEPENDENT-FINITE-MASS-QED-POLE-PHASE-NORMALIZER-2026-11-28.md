# Absolute phase factor for the timelike pole coefficients

The [timelike full-cosine coefficients](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-ANGULAR-2026-11-28.md) retain the squared matrix element, the retarded positive-frequency residue R in ω and the factor 1/V from the charged transfer. They use the angle average dc/2 and leave the k² dk measure and an overall common factor out. At fixed scalar momentum the Boltzmann phase measure for a single charge orientation fixes that factor.

Integrate the spatial delta function over the daughter. Integrating the remaining energy delta over the parent polar angle gives E_b/(pV). The remaining charged phase measure is dE/(8π V) after averaging the parent azimuth. The pole contribution to the photon measure is k² dk dc R/(4π²), after averaging over an isotropic bath and integrating the photon azimuth. Multiplying the external scalar projection 1/(2E_s) gives

```text
k² dk dc dE R / (64 π³ E_s V)
 = [1/(32 π³ E_s)] k² dk (dc/2) dE R/V.
```

Consequently multiply the archived raw-angle coefficients by 1/(32π³ E_s), explicit g² y² for the two bare couplings and k² dk. If the bracket is divided by n_s, undo that division for the unnormalized scalar Boltzmann source. On the source side the derivation refers to a single charge orientation; add a separate charge-conjugate contribution only under the corresponding physical mirroring assumption. The charges already included in the polarization medium are not an additional final-state degeneracy.

In the vacuum specialization R=1/(2k) and ω=k. The expression becomes k/(128π³ E_s V) under dc dE, matching the independent [closed vacuum absorption integral](INDEPENDENT-FINITE-MASS-QED-ABSORPTION-CROSSES-NULL-TRANSFER-2026-11-27.md). The direct gamma-matrix cross-check uses explicitly resolved near-collinear energy neighborhoods in the lightlike limit; a plain Gauss rule across the first occupied bin fails there even as its order rises. For the third archive at k=.00256, an additional geometric angular collar around entering energy-cell edges is needed. Splitting the azimuth alone did not remove the difference. The collar comes from the radial collinear scale divided by the analytic slope of the opening energy with angle; it reduces the checked error to 3.26e−10. Across all four fine and third lightlike high-mesh checks the worst relative difference is 3.61e−6. These extra points do not change the default timelike calculation.

Run `python -W error finite_mass_qed_pole_phase_normalizer_reference.py fine` or replace `fine` with `third`. The tabulated values verify the phase convention for that bare transfer, not the physical replacement of an interacting lightlike cut by a timelike pole.

See the [numerical receipt](results/finite-mass-qed-pole-phase-normalizer-2026-11-28.json).
