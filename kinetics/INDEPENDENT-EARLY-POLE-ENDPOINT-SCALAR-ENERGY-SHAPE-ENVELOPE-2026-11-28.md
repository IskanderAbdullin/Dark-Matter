# Knowing every shell's mean energy does not identify the scalar spectrum

The [energy-conditioned envelope](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR-ENERGY-MOMENT-ENVELOPE-2026-11-28.md) shows that shell counts and mean energies constrain smooth pressure and inverse-energy moments. Test what those same data allow for the cold momentum distribution. As before, allow arbitrary nonnegative measures in each cell, at fixed mass, and normalize either the whole retained count or the cold cells whose upper edge is at most q=1e−5. The cut is an exact shell edge for the four tested charts; in the atomic constructions, boundary atoms retain their assigned cell labels.

In one cell, use the increasing excess energy r=E−m with support [A,B] and mean μ. For its conditional CDF at s, positivity and the first moment give the bounds
```text
if s < μ: 0 ≤ F(s) ≤ (B−μ)/(B−s)
if s > μ: (s−μ)/(s−A) ≤ F(s) ≤ 1 .
```
These come from replacing the energy on each side of the cut by its permitted endpoint. Multiply by cell count fraction and add the lower closed cells. The integral in q of this CDF width is an upper bound on the distance W₁ between any two such measures. Pointwise extremes need not be reached by the same two measures simultaneously. Independently construct a feasible pair: one atom at the momentum corresponding to μ, versus two atoms at the cell endpoints with energy-weighted proportions that reproduce μ. Both have identical shell count and mean energy; their W₁ is the count-weighted sum of the two distances from the mean-energy location to the endpoints.

| Chart | Cold W₁ of explicit same-moment pair | Cold W₁ upper from the CDF width | Cold W₁ upper using count and support only |
| --- | ---: | ---: | ---: |
| Native | 6.817e−7 | 9.440e−7 | 1.376e−6 |
| S2 | 3.324e−7 | 4.606e−7 | 6.663e−7 |
| S4 | 1.663e−7 | 2.305e−7 | 3.328e−7 |
| S8 | 8.317e−8 | 1.153e−7 | 1.664e−7 |

The S8 construction witnesses room for within-cell cold redistribution even when these two moments are prescribed, while the [measured S4→S8 reconstructed step](INDEPENDENT-EARLY-POLE-ENDPOINT-SCALAR8-MATCHED-LADDER-2026-11-28.md) is 8.879e−9. Those are different comparisons: do not add them as an error budget or infer that an actual distribution attains the witness. On full retained support the S8 CDF-width bound is much wider (2.778e−6) because the coarse high-q cells have not been refined. Atoms are mathematical examples, not proposed physical spectra. No prediction outside the archived support or for CLASS follows.

From `kinetics/`, run `early_pole_endpoint_scalar_energy_shape_envelope.py` with native, S2, S4, S8 matched x=.05 states. The [receipt](results/early-pole-endpoint-scalar-energy-shape-envelope-2026-11-28.json) evaluates the CDF bracket integral at two quadrature orders after splitting at the mean-energy threshold. The inequalities are analytic; the quoted integral values are evaluated numerically and are not directed-rounding certificates.
