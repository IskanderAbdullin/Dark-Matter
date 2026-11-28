# One transverse pole occupation cannot balance two archived scalar modes

The [prior angular coefficient test](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-ANGULAR-2026-11-28.md) holds one scalar momentum. Extend its transverse pole evaluation at fixed k=.0015 to eight scalar magnitudes through q=.0008 on both frozen archives. The outgoing scalar energy and momentum are now varied consistently. Its occupation comes from the corresponding archived scalar bin; q=0 is a geometric rest-limit control only, with no extrapolated occupation. The rest-limit coefficient is angle-independent to floating-point accuracy.

At each q the angle-averaged source of that scalar cell, divided by its frozen occupation n_s(q), is still
```text
J_q(n) = n[(1+1/n_s(q)) A_q−B_q] − B_q.
```
Here n is the occupation of the **same** isotropic transverse photon mode at k=.0015, not a parameter that can differ from one scalar momentum to another. The tested slopes are positive at the two endpoints below. Both parent and daughter energy boundary graphs are recomputed as q changes, including when their admissible angular apertures open separately.

| Archive | Balancing n at q=4.1e−6 | Balancing n at q=.0008 |
| --- | ---: | ---: |
| fine | 2.93830e−9 | .1239291 |
| third | 2.39231e−9 | .1169672 |

In both archives there is therefore a common trial interval (2.93830e−9, .1169672) within which this pole gives a positive source for the cold scalar cell and a negative source for the higher-momentum cell. For example n=1e−4 verifies the two opposite signs directly in the raw conditional coefficients. The two balancing points differ by more than forty million in this test; they cannot be replaced by one balancing constant for both cells. This reflects both the changing charged aperture and the measured n_s, not an estimate of the photon population. The actual production integrated over scalar or photon momentum, and the sign of the full matched source, remain unmeasured by this control.

Joint angular/energy/azimuth mesh refinement and the analytic angle-opening checks are retained for every sample. The q=.0008 fixture uses a separately refined energy and cosine mesh; do not reuse the least-expensive cold setting merely because it is accurate at q=4.1e−6. Neither the rest-limit nor the highest-momentum sample is used to infer scalar phase-space weights for an integral.

From `kinetics/`, run:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_timelike_pole_scalar_momentum_reference.py fine > fine-q.jsonl
python -W error finite_mass_qed_timelike_pole_scalar_momentum_reference.py third > third-q.jsonl
python -W error finite_mass_qed_timelike_pole_scalar_momentum_receipt.py fine-q.jsonl third-q.jsonl
```
The [machine receipt](results/finite-mass-qed-transverse-two-scalar-modes-2026-11-28.json) includes independent numerical mesh checks and the conditional trial signs. No population fit, cosmological prediction, or transfer outside the timelike sample is implied.
