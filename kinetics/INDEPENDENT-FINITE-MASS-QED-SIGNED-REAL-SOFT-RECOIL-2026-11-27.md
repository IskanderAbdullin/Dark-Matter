# Exact signed real recoil with a leading thermal soft subtraction

Use the [two real-emission diagrams](INDEPENDENT-EARLY-POLE-YUKAWA-REAL-RADIATIVE-WITNESS-2026-11-27.md)
on both sides of a signed photon energy omega: positive for outgoing
real emission and negative for incoming absorption. This is the
crossed vacuum amplitude. The parent is at rest in the deliberately
chosen photon bath. Integrate the exact two-body recoil of the
daughter–scalar pair, at invariant mass squared M²−2M omega.
The two scalar masses come from the archived pole receipt, but
neither archived occupations nor its plasma boost are assumed.

Factor omega² from the coherent spin sum before evaluation, then
include the recoil phase-space ratio relative to the Born point.
The remaining angular integral G(r), r=omega/M, is smooth at zero.
The independent explicit spin trace agrees with the scaled invariant
form to 4.7e−15 relative on both signs. An angular change of variable
in log(1−beta z) resolves the forward peak. G(0) matches the
[analytic massive radiator](INDEPENDENT-FINITE-MASS-QED-SOFT-SUBTRACTION-2026-11-27.md)
on both starts. The residual integral uses
`G(r)−G(0) + n_B(M r)[G(r)+G(−r)−2G(0)]` divided by r.
For r≤1e−4 use high-precision Taylor coefficients to retain the
small differences rather than subtracting noisy angular integrals.
The direct 65-digit audit of the small-r representation agrees to
1.4e−16 relative.

At r_max=.01 and the specified reference bath T=.9, the vacuum
piece is about −.09250 while the even thermal piece is about
+3.3799 in the normalized integrals before the common gauge factor.
Taking the infrared cutoff from 1e−6 to 1e−10 approaches the
zero-cut result linearly; at 1e−10 the remaining error is about
3.29e−8. The 28/52 angular replays differ below 6.2e−10.
The finite parts depend on r_max and the arbitrary reference bath:
only the leading virtual soft counterterm has been subtracted.
No finite virtual remainder, matching to a hard cut, dressed
widths, or physical thermal decay rate follows.

Replay:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_signed_real_soft_recoil_reference.py
```
[Raw receipt](results/finite-mass-qed-signed-real-soft-recoil-2026-11-27.json).
