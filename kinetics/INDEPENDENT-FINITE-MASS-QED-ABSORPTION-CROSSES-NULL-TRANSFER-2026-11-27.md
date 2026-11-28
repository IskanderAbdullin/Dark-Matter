# A real absorption chart remains finite through the null neutral transfer

For $P+\gamma\leftrightarrow b+s$, hold the same cold scalar momentum
and choose photon energy and angle to that scalar as the outer
variables. The transfer between the two charged lines is
$q-k$, not $q+k$ as in the [emission process](INDEPENDENT-FINITE-MASS-QED-BARE-EMISSION-MONOTONE-BRIDGE-2026-11-27.md).
Its invariant $U=m_s^2-2k(E_s-q\cos\theta)$ crosses zero within the
tested band. An emission-only monotone parent chart therefore
cannot determine the absorption support by a sign extrapolation.

Let $e$ and $V$ be the signed energy and spatial momentum magnitude of the neutral
transfer and $D=M_P^2+U-m_b^2$. The parent aperture comes from

$$
4V^2(E_P^2-M_P^2)-(2eE_P-D)^2\ge0.
$$

For the positive-energy timelike part it gives a bounded interval.
For spacelike transfer the positive-energy support lies above a
finite lower root; the apparent upper branch is not physical.
Rationalizing the lower root avoids cancellation when $U$ reaches
zero. It agrees with the exactly null expression through both
sides within 2.27e−14 relative in the targeted witnesses.
The archived occupation caps are applied separately to both branches.

Average the two charged-photon factors over the exact azimuth,
resolving the narrow charged collinear features in parent energy
and the exact angular knots at archive edges. The azimuthal
comparison with independent quadrature is below 1.19e−11 relative.
The two finer rules for eight sampled photon energies in each
archive agree within 2.78e−5 relative in the positive gain and
inverse weights, before applying the photon Bose factors.
The chart reproduces the independently published signed parent
absorption kernel at .00016 and .0002 to at most 1.45e−5 relative.
Those are shared one-charge real weights, not the signed scalar
collision integral after combining the real and virtual terms.

This chart is a finite-kinematics tool for extending a local
leading-soft subtraction. It is not an inclusive correction,
a resummed photon propagator, or a physical thermal rate.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_absorption_transfer_reference.py
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_absorption_transfer_audit.py
~~~
[Samples](results/finite-mass-qed-absorption-transfer-sampled-2026-11-27.json),
[coordinate gates](results/finite-mass-qed-absorption-transfer-audit-2026-11-27.json).
