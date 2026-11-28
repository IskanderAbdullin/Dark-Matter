# Algebraic audit of the scalar-momentum quadrature knots

The [cold scalar-band quadrature](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-COLD-SCALAR-BAND-2026-11-28.md) enumerates outer scalar q knots using two sign-bracket resolutions. Sign scans alone do not rule out two roots entering a single bracket. This audit solves the finite edge-crossing equations algebraically at the same fixed transverse pole, k=.0015, for both archived distributions.

At an extreme cosine c=±1, all momenta at the charged support boundary are collinear. When the parent histogram edge fixes P, the daughter is `b=P+K−q`; when the daughter edge fixes b, the parent is `P=b−K+q`. Fix either direction for the known charged momentum and set the other charged invariant to its squared mass. In both cases this gives `a sqrt(m_s²+q²)−Bq=C` for known a, B, C. The candidate equation is
```text
(a²−B²)q² − 2BC q + (a m_s−C)(a m_s+C) = 0 .
```
Use both roots or the linear limit, then reject candidates outside the open interval, invalid for the original unsquared equation, below either positive-energy charged shell, or inconsistent with the collinear momenta. The past-null boundary is another linear instance of this equation. Retain the two interval endpoints and exact scalar histogram edges separately. A synthetic round-trip recovers two roots in a bracket with the same sign at both endpoints.

For the cold interval, the following gap and shoulder intervals, and their complete combined interval, every enumerated knot agrees with the corresponding numerical search to 2.43e−18 absolute or better on the replay platform. The combined interval [1e−6,3e−4] has respectively 122 and 170 knots. No additional positive-energy charged orientation enters the interval; each archive has one past-null crossing in the shoulder. The scan agreed here, but its agreement is now checked against the polynomial degree rather than merely another scan spacing.

Run `OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_scalar_q_knot_exact_pilot.py` in `kinetics/` and compare the [receipt](results/finite-mass-qed-scalar-q-knot-exact-2026-11-28.json). This audits support splits within q<k for one equilibrium-free frozen photon pole; it does not certify quadrature error for arbitrary archived changes, cover q≥k, or predict either an actual photon population or a matched rate.
