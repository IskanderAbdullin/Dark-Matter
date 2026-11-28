# Transport the heavy parent without interpolating its phase

The [fixed-background finite-parent check](INDEPENDENT-WARM-LATE-FINITE-PARENT-TURN-2026-11-26.md)
evolves occupations at native quadrature nodes. Those physical nodes
change with the daughter momentum support under expansion. Copying
occupation by array index would silently change the represented
phase.

Build a common comoving parent grid from the union of the daughter
geometry breakpoints at selected scale factors. Keep the same
positive number weights and parent occupations during collisionless
expansion. At each collision epoch, reconstruct the lab energy
E = sqrt(M² + Q²/a²) from the fixed Q coordinate and re-cut the
daughter angular panels at that epoch. This retains the original
angular and mass kinematics, restricted to the selected closed
daughter support. The parent momentum used in the decay weight is
computed from Q/a, rather than recovering a nearly zero momentum by
subtracting the rest energy from E.

For the .68 and .80 GeV tested late phases, at the initial, middle
and final scale factors of their short local windows, the order-4
common chart matches the native signed W current to 1.12e-10
relative or better. At order 6 the difference is below 2e-14.
The stationary integrated number and energy invariants hold to
better than 3e-17 relative to the selected gross number rate in
these units. This validates geometry on a prescribed phase, not
the extrapolation of the original fixed-background trajectories.

In the collisionless parent drift comparison, compute pressure-trace
work by integrating d(a⁴ρ)/dlog(a). An independent stable difference
of sqrt(a² M² + Q²) at the endpoints verifies this work without
subtracting two nearly equal physical energies. Initial work for the
selected .68 short window is 9.46e-28 GeV⁴ and for the .80 window is
5.68e-26 GeV⁴. Number and the shape at fixed Q are invariant by design.
The separate a³ρ redshift energy loss is independently checked using
the integral of −p²/E and its cancellation-free endpoint difference.
As a stress test, extend the prescribed .80 drift to Δlog a=.05.
The signed stationary current on the union chart still agrees with
the native chart within 2.40e-6 relative across three sampled scales;
the pure-parent a³ energy loss is −6.720e-22 GeV⁴. Neither control
refits a parent distribution or assumes an exterior return channel.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MALLOC_ARENA_MAX=2
# Prepare turn6 and competing080/tapes from the corresponding notes.
for n in 4 6; do
 python -W error warm_comoving_parent_chart_audit_reference.py EOS.dat turn6/step-4.json "chart068-$n.json" --order "$n" --span 1.5e-8
 python -W error warm_comoving_parent_chart_audit_reference.py EOS.dat competing080/tapes/step-5.json "chart080-$n.json" --order "$n"
done
python -W error warm_comoving_parent_chart_audit_reference.py EOS.dat competing080/tapes/step-5.json chart080-broad.json --order 6 --span .05
~~~

[Raw geometry and collisionless-work diagnostics](results/warm-comoving-parent-chart-audit-2026-11-26.json).
