# A common plasma soft cut in two photon angle charts

Keep the [exact crossed and eikonal-subtracted real recoil](INDEPENDENT-FINITE-MASS-QED-SIGNED-REAL-SOFT-RECOIL-2026-11-27.md)
but put the photon occupation and the unresolved photon energy cap
in the **plasma** frame. Integrate all scalar angles, so this does not
evaluate the cold scalar-momentum bins. The two parent momenta
(archive kinematic midpoints) correspond to boosts gamma ≈28.57
and 28.73. Use the same illustrative plasma temperature T=0.9.

If xi is the rest-frame photon angle to the boost, the lab energy
is D omega_rest with D=gamma(1+beta xi). For the same physical
lab cutoff and occupation, integrate both over isotropic rest
angles and over plasma photon angles with their aberration
Jacobian. Use the exact recoil matrix integral, with the prior
Taylor coefficients only for the cancellation-sensitive soft
differences. The two angle/quadrature representations agree
within 8.5e−14 absolute across the two poles and the
2.5e−6, 5e−6, and 1e−5 plasma photon caps.

There is a useful small-cap check. The leading vacuum recoil
is proportional to the rest-angle average of D⁻¹, about
0.1417 and 0.1411 on the two starts. The leading even thermal
recoil is proportional to the corresponding D⁻² average, which
is exactly one. At the 5e−6 cap the vacuum residual is about
−0.000132 and the even thermal residual about +0.033804 in
the normalized integrals before the gauge factor. This is a
restricted all-angle limit, not smallness of the actual
cold-bin result. The finite one-loop virtual terms, dressed
propagators, quantum statistical occupations of the fermions,
and cutoff matching to the hard real operator remain open.

Replay from kinetics with the earlier soft and signed-recoil receipts:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_boosted_soft_recoil_reference.py
```
[Raw receipt](results/finite-mass-qed-boosted-soft-recoil-2026-11-27.json).
