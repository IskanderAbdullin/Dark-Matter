# Single-cut causal assignment for the scalar-latest triangle

The [bare absorptive reference](INDEPENDENT-FINITE-MASS-QED-CHARGED-ABSORPTIVE-GUARD-2026-11-27.md)
used the same positive imaginary guard for each charged denominator.
Here fix a definite retarded three-point ordering before comparing
the two positive charged cuts. This follows the branch-sum
construction of [Brandt et al.](https://arxiv.org/abs/hep-th/9809177).

Orient three internal edges consecutively around the loop. Fix
one vertex to the plus branch as the latest operator, and sum the
two other branch labels including the two vertex signs. For each
internal propagator write the common thermal addition as T_i
independent of its branch labels. A symbolic expansion cancels
all products of two or three different T_i. With latest vertex 0,
the three surviving thermal coefficients are

    T_0 A_1 A_2,     T_1 R_0 A_2,     T_2 R_0 R_1.

All three choices of latest vertex give the corresponding rotated
assignment. R and A denote vacuum retarded and advanced denominators;
this is a branch algebra statement, not a normalization of a
collision rate.

Choose vertex 0 as the scalar insertion, with oriented lines
0,1,2 carrying daughter b−K, photon −K, and parent P−K.
In the positive occupied radial overlap, the remaining daughter
energy in the parent cut, and the remaining parent energy in the
daughter cut, are both positive. They therefore use +i0 (R) and
−i0 (A), respectively. These imaginary signs oppose before any
common charged statistical factor is applied. The separate real
PV parts do not change sign from this assignment. This statement
only covers the positive near-coincidence charged cut pair, where
the uncut photon is spacelike at nonzero displacement.

For the charged leg with on-shell bath-frame energy E, momentum p,
mass m, and cosine c relative to the scalar momentum, the local
width and root slopes A,B satisfy

    Q p² sqrt(A²+B²) = sqrt[(P_leg·Q_scalar)²−m² m_scalar²].

The right side is identical for the parent and daughter on the
same Born event. Hence the scalar-latest imaginary log coefficient
in this overlap, with no common overall statistical sign or
coupling attached, is proportional to f_parent−f_daughter, with
one shared invariant coefficient. This identity holds on each
continuous occupied cell, not at an unspecified bin assignment.
The exactly opposite local pieces would cancel for equal occupations.

At three interior points per archived grid the parent/daughter
local slopes are +.448623/−.128054 (fine) and +.528075/−.285476
(third), leaving +.320569 and +.242599 in the chosen convention.
The actual finite-guard radial slopes agree to 1.4e−7 absolute.
This comparison is not a sign or magnitude of physical damping:
the photon cut, other branches, and matched real/external response
must still be organized for the physical observable.

Replay with SymPy, NumPy and SciPy from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_single_cut_contour_reference.py
```
[Raw numerical receipt](results/finite-mass-qed-single-cut-contour-2026-11-27.json).
