# Timelike pole coefficients over all photon angles at fixed k

The preceding [fixed-angle calculation](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-AFFINE-SOURCE-2026-11-28.md) leaves the mode occupations unspecified. Extend the same bare retarded-residue and Pauli coefficients to the full cosine range at fixed scalar q and fixed pole k. Both endpoint histograms are held fixed. The normalized angle measure is dc/2; the remaining k² dk and common couplings and phase constants are not included.

## Exact opening between histogram cells

At transfer t=q⁰−ω<0, a parent at energy E has momentum p and the daughter at E−t has momentum b. For transfer magnitude V the angular root exists precisely when |b−p| ≤ V ≤ b+p. Converting V²=k²+q²−2kqc gives an exact cosine at each entering parent or daughter histogram edge. In this soft fixed-k sample only the |b−p| edge is relevant. The implementation checks that no upper edge was discarded and splits the possible past-null boundary too. The independently solved radial aperture agrees at the collected angle knots to 1.19e−14 relative.

The batched radial implementation uses the on-shell numerator reduction
```text
(2 P·ε + slash(K) slash(ε)) / ((P+K)²−M²)
+ (2 b·ε − slash(ε) slash(K)) / ((b−K)²−m²).
```
Spinors apply the same two timelike polarization conventions and their physical retarded residues. A separate uncontracted spinor expression with direct denominators checks the batched matrix result against the full numerator to 6.29e−9 on the tested angular points, and the retained independently implemented radial quadrature to 1.29e−8. Full angular integration jointly refines the cosine quadrature, energy and azimuth; maximum coefficient mesh difference is 2.86e−7 relative.

## Open forward and inverse support separately

At k=.00135 the transverse mode gives zero on both archives for this scalar momentum. Fine transverse k=.001395 and .001405, and third k=.001415 and .001425, have a nonzero full-angle forward coefficient but exactly zero inverse support. For example fine at .001405 has A≈1.3350e−15 and B=0 in the common omitted units. Discarding angle bins where either charged archive is empty would incorrectly remove these genuine forward contributions. Their smallness also means that being open does not make them dynamically dominant.

At transverse k=.0015 the cosine average increases A by 27.46% (fine) and 39.91% (third) compared with the cosine-zero value. This is near the exponentially sparse tail opening. The balancing pole occupation for that full-angle average is then 2.9383e−9 and 2.3923e−9 respectively, assuming a common occupation over directions. The leading relation is still n[(1+1/n_s)A−B]−B when n at this k is angle-independent. A non-isotropic pole population would have to be retained inside the angle integral. This fixed-k result supplies neither the actual n nor a momentum-integrated scalar source.

From `kinetics/`, run the following three references, whose datasets form the [combined receipt](results/finite-mass-qed-timelike-pole-angular-2026-11-28.json):
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error finite_mass_qed_timelike_pole_batched_absorption_reference.py
python -W error finite_mass_qed_timelike_pole_angular_reference.py fine
python -W error finite_mass_qed_timelike_pole_angular_reference.py third
```
