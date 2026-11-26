# Evolution to 12 MeV and continuum common-state capacity

The 50 MeV benchmark is continued to the scale factor corresponding to
12 MeV on the reference visible adiabat. Visible feedback and the exterior
ledger are retained, and the first dark momentum cell starts at zero.
At 48 steps the final scalar energy share is 0.2268980 with vacuum source
and 0.2268993 with the quantum source, versus 0.2272766 at 24 steps with
vacuum. The quantum run reaches actual visible temperature 11.99971 MeV.
It gives represented scalar and stable yields 6.14329e-6 and 1.59772e-5
relative to initial visible comoving entropy. The scalar number fraction
below q/T_R=0.01 is 1.53e-7. The represented dark energy and count ledgers
close within 9.1e-13. These are numerical finite-time controls, not relic
abundances and not a resolved singular state at zero momentum.

For a separate common equilibrium diagnostic at the measured represented
counts, solve for T_c where continuum scalar excited states at mu_phi=m_phi
carry the full scalar count. At T_c match the stable fermion chemical
potential to its own count, then compute the total kinetic energy above
the two rest contributions. The fixed-number kinetic energy increases
with common T: each derivative is the kinetic variance after removing
the count fluctuation. At lower kinetic energy the ideal common state
would allocate part of the scalar number to the ground mode. This
requires isolated common equilibrium, not merely elastic collisions.

The controls verify the massless saturation formula, nonrelativistic
asymptote, and manufactured normal and condensed examples. Each 45, 90,
20, 15 and 12 MeV snapshot allows a normal common state. At 12 MeV the
kinetic energy is roughly 19 times its saturation value and the scalar
fugacity relative to the rest threshold is about 3.6e-4. The actual
retained shell entropy and spectra are separately compared with the
common continuum state without imposing it on the kinetic solution.

From `kinetics`, with the external central EOS:

```sh
export OPENBLAS_NUM_THREADS=1
for SRC in matched-vacuum-order quantum-visible-source; do
 python -W error joint_resonant_frw_split_reference.py EOS.dat --physical-amplitude --visible-source --"$SRC" --visible-backreaction --trace-exterior --with-elastic --empty-pulse 1e-8 --empty-stop .24 --support-low .005 --steps 48 --vacuum-scaled-residual --allow-empty-shells --include-zero-cell --save-snapshots snapshot-"$SRC".json
done
python -W error continuum_dark_normal_capacity_reference.py EOS.dat results/12-quantum-48-2026-11-26.json --snapshots results/12-quantum-48-shape-2026-11-26.json
```
