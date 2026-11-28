# Two coalescing charged poles at an archived occupation edge

The [separated occupation-edge calculation](INDEPENDENT-FINITE-MASS-QED-HISTOGRAM-EDGE-POLE-2026-11-28.md) has only one sequential pole inside the diagnosed window. The [earlier merging-pole control](INDEPENDENT-FINITE-MASS-QED-CUT-POLE-COALESCENCE-2026-11-27.md) kept coherent residue cancellation but did not move a step into the two-pole region. Here vary both locations together in actual archived spacelike shell coordinates. No physical width sign is chosen.

Let x=D1/ε, h=2 q·K/ε, and put the occupation step at x=η. At equal real residues the squared coherent kernel with the same imaginary orientation is
```text
F(x;h) = (h²+4) / [(x²+1)((h−x)²+1)].
```
Its whole-line integral is 2π independent of separation. The leading left cumulative primitive is
```text
I(η;h) = π+atan(η)+atan(η−h)+log[(1+η²)/(1+(η−h)²)]/h.
```
The log quotient at h=0 is replaced by 2η/(1+η²). For a step with one-sided weights w⁻,w⁺ the scaled leading area becomes w⁻ I + w⁺ (2π−I). Reversing the second imaginary orientation multiplies this leading expression by h²/(h²+4). These identities follow from the rational partial fraction and do not assume a cancellation between positive probabilities.

The numerical replay changes ω so that the measured h is held fixed while ε is varied; it separately solves the azimuth so an actual occupation edge lies at the requested η. There are 60 tests on the two archived spectra, with 5 separations, 3 edge positions and 2 small widths at k=.0008. A separate numerical integral verifies the dimensionless cumulative formula to 8.89e−16 absolute. Full spinors with actual occupation weights and physical shell Jacobian reproduce the asymptotic prediction to 1.89e−5 absolute; their relative residue-difference squared is below 7.80e−12.

At centered edge η=0 and h=2, the third archive's common-sign leading area is about 7.17910, not 2π, and the opposite-sign case gives about 3.58955. At h=0 the leading opposite-sign area vanishes; finite normalized spinor remnants remain near 2.4e−11 in the observed fixtures. Neither removing the step nor multiplying a smooth result by a constant edge value recovers the general joint limit.

Run `python -W error finite_mass_qed_coalescing_histogram_edge_reference.py` from `kinetics/` and compare the [receipt](results/finite-mass-qed-coalescing-histogram-edge-2026-11-28.json). This is a mathematical diagnostic of a regulated bare overlap on real archive geometry. It is not a rule for a medium width, a physical sequential subtraction or a matched scalar rate.
