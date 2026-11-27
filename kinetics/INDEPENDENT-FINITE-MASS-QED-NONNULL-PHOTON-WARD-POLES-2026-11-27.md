# A spacelike photon cut opens genuine on-shell charged intermediate states

In the [frozen transverse Landau cut](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-CUT-2026-11-27.md),
test the two bare charged diagrams for $P+K\leftrightarrow b+q$
with a non-null photon momentum $K$. Preserve the exact two
inverse propagator denominators $D_1=K^2+2P\cdot K$ and
$D_2=K^2-2b\cdot K$. Their sum is exactly $2q\cdot K$.
The null-photon propagator shortcut would lose this relation.
Use explicit positive spinor sums with two spatial transverse
polarizations. Shifting either by a multiple of $K$ changes the
sampled off-pole result by at most 2.47e−9 relative. This is a
bare Ward consistency test, not a finite-width gauge prescription.

For the sampled timelike transverse pole, neither charged
denominator crosses zero. For each of the sampled spacelike
geometries and azimuths both do. As an example at momentum
.0008 and frequency .0002, the fine case with both specified
angles zero reaches poles at parent energy 1.42198497 and
1.42218380. In the independent reconstruction all 48 roots
have future-directed internal charges and lie in occupied and
unblocked external archive cells, with maximum absolute internal
mass-shell error 1.15e−15. Thus a pole is not an artifact of
continuing into an empty region of the frozen spectra.

These are **sequential kinematic overlaps**. One cannot multiply
the retarded spacelike spectral weight by the raw squared bare
amplitude and interpret the result as a new finite collision
channel. Consistent treatment of the reducible part, incoming
statistical photon correlator, and charged/photon dressing are
still required. In particular no physical width is inferred from
the distance between two kinematic roots.

Replay in kinetics with the pinned archived spectra:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_nonnull_photon_ward_poles_reference.py
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_nonnull_photon_occupied_roots_audit.py
~~~
[Ward and pole locations](results/finite-mass-qed-nonnull-photon-ward-poles-2026-11-27.json),
[future/occupation check](results/finite-mass-qed-nonnull-photon-occupied-roots-audit-2026-11-27.json).
