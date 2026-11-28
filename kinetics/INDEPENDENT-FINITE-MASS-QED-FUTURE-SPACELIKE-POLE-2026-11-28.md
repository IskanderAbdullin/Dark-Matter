# Continue the frozen transverse pole calculation through positive transfer energy

The [finite-q support audit](INDEPENDENT-FINITE-MASS-QED-FUTURE-Q-SUPPORT-2026-11-28.md) shows that a timelike photon momentum K does not imply a timelike charged transfer T=q−K. The previous radial routine only accepted T⁰<0. For the k=.0015 transverse fixtures, the pole is about twice as massive as the scalar, so the only potentially allowed charged transfers are spacelike on either side of T⁰=0. This change explicitly opts into the additional branch; default past-directed calls retain their existing behavior.

Let t=T⁰, v=|T spatial|, s=T²<0, and A=M²+s−m². The lower parent-energy limit which enforces the real transverse momentum is
```text
E_min = (v sqrt(A²−4M²s) − A t) / (−2s).
```
For A t>0 use the algebraically equivalent `(A²+4M²v²)/(2(v sqrt(A²−4M²s)+A t))` to avoid near-null cancellation. The test reconstructs the same limit by boosting from the zero-energy frame of T at 65 decimal digits. Across 15 fixtures including near-null separations of 3e−13, the maximum relative defect is 5.78e−16. On the positive-energy branch the angular partition also permits the collinear boundary where the two charged spatial magnitudes add rather than subtract.

In the two archived distributions the transfer energy changes sign around q=.0021328 and .0021359. At fixed cosine zero the three narrow brackets show that the affine coefficients converge continuously across this crossing: offsets ±1e−7 give a maximum one-sided contrast of 3.47e−4, reduced to 3.47e−8 at ±1e−11. The opt-in past-side values equal the existing path exactly in those fixtures. The spinor trace is checked against the unreduced vertex on both sides.

For q = .004, .015, .03 and .08, integrate the full cosine at fixed k and check each of the independent affine coefficients across a mesh ladder. The worst accepted componentwise difference is 1.93e−5. The q=.03 fixtures exercise the additional angular boundary. The coarse q=.08 results did not pass and were refined. The conditional scalar-number-neutral n at the first q is .00156568 / .00165302 (fine / third), and at q=.08 it is 6.35748e−7 / 7.00587e−7. They are much larger than the previously resolved cold-band neutrality near 2–3e−9. This is a statement about sampled q at one fixed k, not about the actual n or the magnitude, sign or neutral point of a momentum-integrated source.

Run `OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_future_spacelike_pole_pilot.py` from `kinetics/` and compare the [receipt](results/finite-mass-qed-future-spacelike-pole-2026-11-28.json). The replay uses the same unbroadened transverse pole, frozen distribution and unresummed charged vertex as the past-side reference. It does not add the Landau continuum, self-energy/vertex dressing, or the other matched collision channels.
