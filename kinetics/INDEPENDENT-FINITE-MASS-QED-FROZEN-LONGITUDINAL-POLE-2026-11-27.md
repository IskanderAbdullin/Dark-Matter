# A longitudinal pole opens the same charged transfer earlier than the transverse pole

In the same [frozen collisionless charged response](INDEPENDENT-FINITE-MASS-QED-FROZEN-TRANSVERSE-PLASMON-GATE-2026-11-27.md),
solve the longitudinal Coulomb-dispersion condition with nonzero
charged masses. Use the same two spin states and mirrored charged
populations. For timelike frequency the longitudinal kernel is
$\operatorname{atanh}(kv/\omega)/(kv/\omega)-1$, weighted by
the charged radial derivative. The root sets the integrated
kernel equal to $k^2$. Direct charged jumps and independent
positive radial integration agree within 4e−13 relative,
including at $\omega=k$.

The longitudinal root has the same $k=0$ plasma frequency as
the transverse root. Its small-$k^2$ slope is .59957 in both
fixtures, consistent with the near-ultrarelativistic value
3/5. With finite constituent masses the maximum particle speed
is below one on the archived domain. The timelike longitudinal
curve reaches the light cone at finite momentum about
.005761–.005764; this is not an extrapolation of an undamped
pole into the spacelike cut.

The longitudinal version of $q-K$ is past timelike at small
momentum and becomes spacelike earlier than for the transverse
pole. On the sampled directions this kinematic opening lies
between .000948 and .000969, with energy around .001824–.001832.
These are neutral-transfer openings, not charged occupation
activation or a production threshold in a complete rate.
For the equilibrium form of the dispersion see
[Braaten and Segel](https://arxiv.org/pdf/hep-ph/9302213).
No longitudinal external matrix element, width or finite
scalar production/absorption rate is computed here.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_frozen_longitudinal_pole_reference.py
~~~
[Raw receipt](results/finite-mass-qed-frozen-longitudinal-pole-2026-11-27.json).
