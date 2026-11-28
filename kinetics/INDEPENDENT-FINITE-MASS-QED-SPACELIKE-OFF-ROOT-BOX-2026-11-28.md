# Finite spacelike scalar and photon-energy exchange without on-shell charged roots

Extend the [bounded fixed-geometry off-root test](INDEPENDENT-FINITE-MASS-QED-SPACELIKE-OFF-ROOT-PATCH-2026-11-28.md) to a finite external box:
`k ∈ [.001265,.001295]`, `ω/k ∈ [.493,.507]`, and `cos(q,K) ∈ [−.54,−.46]`.
The same two parent energy intervals, `[1.18,1.26]` and `[1.85,2.25]`, are kept separate throughout. This is not a slicing prescription for a total collision integral.

Outward-rounded interval operations certify spacelike charged transfer, a real parent aperture, and nonzero charged denominators over every point of each closed box and every parent azimuth. The latter enters the interval only through cosφ in [−1,1]. The smallest certified denominator magnitude is 1.13683e−4, including variation of k, ω/k, geometry and parent energy. No inverse width, excised resonant boundary integral, or interleaved cut subtraction is invoked.

The [one-charge phase normalization](INDEPENDENT-FINITE-MASS-QED-POLE-PHASE-NORMALIZER-2026-11-28.md) now includes the continuous dω/(2π) in place of the isolated pole residue. In u=ω/k coordinates the integration uses k³ dk du, the cosine average dc/2 and the fixed scalar factor 1/(32π³E_s). The positive-frequency leading continuum commutator and its charged noise set the directional coefficient; no free photon Bose population is added. Quadrature separately splits every archived charged velocity edge inside the u box, including the longitudinal logarithmic feature. The checked maximum mesh change of the integrated coefficients is 3.27e−8 relative.

Summing the two parent patches, the following coefficients retain their physical phase weights but are divided by n_s and omit explicit g²y². The scalar q is held at 4.1e−6. Photon heat is the signed photon energy gained or lost, so the net energy carried from this photon sector is recorded with the opposite sign from the scalar-number event.

| Frozen archive and mode | Scalar number coefficient | Photon energy coefficient |
| --- | ---: | ---: |
| fine transverse | +1.805953e−8 | −1.155877e−11 |
| fine longitudinal | +7.968721e−9 | −5.101243e−12 |
| third transverse | +1.910814e−8 | −1.222992e−11 |
| third longitudinal | +8.413454e−9 | −5.385945e−12 |

A separate shape-only sensitivity toggles the previously computed exact free finite recoil in the statistical ratio on the original fixed-geometry patches, keeping the retarded commutator fixed. The positive signed coefficients change by at most 6.96e−6 relatively and reconstruct the same gain-minus-loss expression to 1.78e−15. This diagnostic does not substitute a finite-recoil retarded propagator or determine a dispersive real part.

From `kinetics/`, reproduce the [box receipt](results/finite-mass-qed-spacelike-off-root-box-2026-11-28.json) and [finite-recoil sensitivity](results/finite-mass-qed-spacelike-off-root-recoil-2026-11-28.json):
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_spacelike_off_root_box_reference.py
python -W error finite_mass_qed_spacelike_off_root_recoil_sensitivity.py
```
The signs are local to these certified boxes. No inclusive spacelike contribution, quantum-width prescription, population evolution, or matched physical scalar production is inferred.
