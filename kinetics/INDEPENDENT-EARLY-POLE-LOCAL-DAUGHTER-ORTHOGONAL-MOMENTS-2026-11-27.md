# Single-moment daughter controls on the complement of the initial aperture

Leave the two prepared daughter aperture cells exact, then use the
same two-parameter affine FD changes in the other 18 cells. Match
one full retained moment to native while holding the other at its
unmatched prepared value. This separates *target constraints* for
one retuning recipe; it does not keep interaction energy or higher
moments fixed.

Matching number while preserving the unmatched bare energy reaches
both targets within the original tolerance. Its affine coefficients
are about +0.00620 and -0.19123, under the magnitude gate 0.2.
The initial source shift changes by about 105 ppm relative to the
unmatched aperture-local shift. The bare-energy-only control reaches
the moment target numerically but needs a coefficient about +0.2513,
so it is rejected by the unchanged amplitude gate and no rate is
reported for it. The paired number-and-bare-energy control still
requires only +0.00616 and +0.05337 and shifts the initial result by
about 33 ppm.

These checks are static. They do not identify a unique decomposition
of the finite compensated effect. Use the paired fine archive with
`early_pole_daughter_local_orthogonal_moment_reference.py` to replay
the proposed constructions and amplitude check.

## Finite held-response control at the same starting number

The accepted static number-only preparation was evolved using the
same unsplit graph, support, initial cold seed and numerical
tolerances as the raw and two-moment controls. The static and
evolved retunes match exactly and the four terminal ledgers pass.

| Aperture-local initial preparation | Finite mean shift from native |
|---|---:|
| Raw substitution | +5.17245309e-4 |
| Number matched, bare energy held at raw value | +1.39448067e-3 |
| Both number and bare energy matched | +8.85152056e-5 |

At fixed native daughter count the two manufactured compensated
inputs have finite shifts differing by 1.30597e-3, even though
their initial cold-rate suppressions are nearly identical. The
number-only shift is about 2.70 times the raw one. This establishes
a sensitivity to the compensation choice, not a unique derivative
with respect to the bare energy: the retunes alter the rest of the
daughter distribution as well, and the energy-only arm is not
available under the declared amplitude gate.

Run the simultaneous chart with the paired archive and
`--number-only-local-retune`; the static and evolved results are
cross-checked by `early_pole_daughter_number_only_time_compare.py`.
