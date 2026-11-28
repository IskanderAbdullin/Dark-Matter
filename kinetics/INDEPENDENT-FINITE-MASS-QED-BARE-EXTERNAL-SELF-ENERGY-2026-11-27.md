# Resolve the bare charged and null-photon external self-energies

After the [longitudinal Ward check](INDEPENDENT-FINITE-MASS-QED-LONGITUDINAL-WARD-2026-11-27.md),
evaluate the separate Feynman-gauge bare on-shell self-energy projections
and their fixed-spatial-momentum, fixed-external-spinor energy derivatives.
The two charged masses are nonzero. These are scalar projections of
matrices with the spin average normalized by 4E. They are not in-medium
spinor projectors or fully matched residues.

For a charged internal four-vector L with energy sign ± and signed
energy displacement Δ=L⁰−E, the angular numerator from the gamma-matrix
contraction is (K²+2m²)/E with K=P−L on shell. For the two elementary
angular integrals J1=∫du/K² and J2=∫du/(K²)² the projected derivative is
2Δ(J1+2m²J2)/E. Its denominator endpoints obey

    K²(u=+1) K²(u=−1) = 4m² Δ².

With radial measure r/2 this reduces the singular piece exactly to
r/(E Δ), not a cancellation between two large numerical values.
For the positive cut the remaining subtraction
(r−p)/(E Δ) is (t+E)/(E(r+p)); the leftover Δ J1 term vanishes at
coincidence. For the negative cut Δ never crosses zero. Gamma-matrix,
resolved angular, and complex-step derivative checks agree to 2.1e−10
relative or better.

On the central events the frozen positive-cut derivative is about
+0.303/+0.417 (fine/third parent) and −0.260/−0.352 (daughter); the
negative-cut values are respectively −2.609/−2.597 and −0.537/−0.525.
The two energy signs are separate, not yet combined into a physical
thermal operator. At selected archive step centers the bare radial
log coefficient is p/E times the occupation jump, with either sign
possible. The actual cutoff slopes reproduce it, and the subtracted
cutoff ladder agrees with the analytic cell treatment to about 3e−11
absolute. At the exact jump this is again a scale-defined finite part
at reference radius 3e−4, not a physical resolution prescription.

For the bare **null photon** cut, pairing K and −K at fixed k and
integrating over the full polar cosine yields

    on-shell projection:    −4/E
    fixed-spin derivative:  −4/k² + 2 log[(E+p)/(E−p)]/(E p).

The finite second term has a 4/E² rest limit. Independent angular
integration and complex-step differentiation verify the formula.
With unsigned Bose weight and radial measure k/2, this finite term
is proportional to the Planck moment π²T²/6. The singular −4/k² term
has explicitly been left out of that quoted finite moment. This is
not a prescription for dropping it in a matched medium rate.

For physical context see the [real-time massive self-energy calculation](https://doi.org/10.1103/PhysRevD.31.3280).
This audit does not assign its real-time statistical signs to the
nonthermal archive, dress the photon, shift the external Born aperture,
or compose the vertex, residues and real rate.

Replay from kinetics with the existing spectra:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_bare_external_self_energy_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-bare-external-self-energy-2026-11-27.json).
