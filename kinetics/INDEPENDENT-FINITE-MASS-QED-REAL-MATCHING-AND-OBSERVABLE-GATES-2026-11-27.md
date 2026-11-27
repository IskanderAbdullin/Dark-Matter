# Check real-band matching and which photon observable is being computed

## Separate charts for the same finite absorption band

Compare the inverse two-body cumulative chart against the signed
parent aperture for the same lab photon energy band [4e−6, 1e−5].
This isolates *real* absorption and its inverse. It does not include
emission or a virtual subtraction. Compare the positive forward
and inverse weights before forming the much smaller matched result.

| Archive | Parent-chart real gain | Parent-chart real inverse loss | Cumulative-chart relative gain difference | Cumulative-chart relative loss difference |
|---|---:|---:|---:|---:|
| fine | 791362.12 | 226403.73 | −9.97e−6 | +2.11e−4 |
| third | 855202.58 | 370190.80 | −1.03e−4 | +1.18e−3 |

The signed-parent photon integration changes by about 1.1e−9 under
the two tested outer rules. Its numbers are raw, not scalar-number
corrections. The cumulative inverse check reproduces their scale
but its loss mismatch is **too large for subtracting independently
rounded raw pieces to infer the finite matched source**. Use common
local real/virtual subtraction instead of splicing the cumulative
real rate directly into the small residual.

## Conditional photon spectrum is not the same as kinetic heat

For a single photon mode, use a canonical geometric initial
occupation N. For each N, give the selected event a positive
emission intensity a(N+1), absorption intensity bN, and an
unshifted intensity c+dN. The actual unconditioned kinetic moment
is the mean signed photon jump, weighted by its event intensity.
If instead we compare the post-selected photon spectrum with the
*unselected* bath, the result has an extra contribution:
the covariance of N and the selected-event intensity divided by
the mean event intensity. The covariance cancels if the comparison
uses the same selected events on both sides.

At detailed balance the adjacent N→N+1 and N+1→N probabilities
match exactly. The actual photon heat then vanishes, although
post-selected minus unselected bath energy need not vanish.
At k=4e−5 and T=1, a positive manufactured fixture has essentially
zero kinetic heat but an apparent shift about 0.99995 in the latter
comparison. Direct discrete oscillator sums reproduce the analytical
identities within 5.1e−13 relative or scaled absolute. This is an
observable-definition check, not a prediction for a plasma spectrum.

[Landshoff and Taylor](https://arxiv.org/pdf/hep-ph/9402315) compute
a spectrum per selected decay in a different process and identify
a term associated with bath fluctuations. The one-mode algebra
here explains why conditioning can produce such an additional term;
it does not supply that paper's complete field-theory expression
for this model, and it is not an extra heat term to add to the
unconditional kinetic collision integral.

Replay from kinetics with the pinned archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_photo_soft_overlap_reference.py
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_photon_selection_covariance_reference.py
~~~
[Finite-band overlap](results/finite-mass-qed-photo-soft-overlap-2026-11-27.json),
[conditional-observable test](results/finite-mass-qed-photon-selection-covariance-2026-11-27.json).
