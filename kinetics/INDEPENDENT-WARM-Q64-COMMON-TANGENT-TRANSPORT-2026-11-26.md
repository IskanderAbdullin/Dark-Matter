# Transport a prescribed common-equilibrium pair tangent through the saved q64 clock

At five saved q64 checkpoints between 8.133100 and 7.394049 MeV,
repeat the hypothetical
[common FD/BE equilibrium fit](INDEPENDENT-WARM-PAIR-COMMON-MOMENT-EQUILIBRIUM-2026-11-26.md)
to the represented aggregate number and energy. Keep the old
collisionless cohorts in H, not in the fitted reservoir. At
each anchor the local pair tangent along the two transfer modes
gives a slow rate 1.050H → 0.812H and a fast rate 6.072H → 5.033H.
All fits are hypothetical local comparisons, not measured
thermal shapes or an assertion that their trajectories exist.

In each local heat basis the restoring matrix is the symmetric
Onsager response times the sum of the two inverse susceptibility
matrices. Convert it back to the same fixed moment variables
(dN, dE / 7.5 MeV) at every anchor. Interpolate the entries of
that **fixed-coordinate** matrix divided by the actual anchor H,
linearly in ln(a). Evolve its homogeneous collision-only tangent
by ordered midpoint exponentials from a=1 to 1.10.
This comparator omits equilibrium driving, pressure work,
dark/visible clock changes under perturbations and all other
collision channels.

The fine ordered propagator has eigen-amplitude factors 0.58950
and 0.91558. Thus even this internally thermal pair-only tangent
is not close to instantaneous relaxation on this short clock.
The two positive eigenvalues at each time also agree under
independent symmetric (matrix root and Cholesky) similarities;
this avoids treating the eigensolve of a formed nonsymmetric
product as the exact reference at machine precision. The
ordered propagator's log determinant agrees with the integral
of the negative trace to within 2.1e-15. The fixed-coordinate
Euclidean norm can grow along a mixed transient direction;
it is not the thermal entropy norm and is not an instability
inferred for the evolving archive.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
# ladder/wide is the work directory from the matched q64 run.
python -W error warm_corrected_common_thermal_anchor_tape_reference.py EOS.dat \
  ladder/wide thermal-anchors >thermal-anchors.json
python -W error warm_common_pair_collision_only_propagator_reference.py EOS.dat \
  thermal-anchors.json ladder/wide
~~~

[Five fitted tangents](results/warm-q64-collision-only-common-tangent-inputs-2026-11-26.json.gz);
[ordered propagator and controls](results/warm-q64-collision-only-common-tangent-2026-11-26.json).
