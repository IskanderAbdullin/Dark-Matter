# A reversible scalar elastic reference with number, energy and entropy checks

Build a positive weak-form collision graph from ordered-pair phase
space, the isotropic identical-scalar Born cross section, and the
two scattering angles relevant to outgoing momentum magnitudes.
The signed weak form has factor one quarter over ordered incoming
pairs: one half counts identical initial pairs and one half allows
the difference between forward and reverse in the same integral.
Project each outgoing kinetic energy onto adjacent nodes with
nonnegative barycentric weights. Every retained collision conserves
two particles and projected energy. Out-of-grid events are removed
from both directions, not interpreted as physical disappearance.

The occupation factor uses a cancellation-safe Bose affinity,
with one-way population updates at exactly empty participating
cells. If S is the signed per-event change vector and η the
discrete entropy gradient, the rate has the sign of S·η.
Then R(S·η) is nonnegative event by event. Thermal detailed balance
at fixed positive mass is exact up to roundoff since η is affine
in energy at a Bose distribution. Fractional population weights
are a numerical closure, not fractional physical particles.

A smooth manufactured trajectory passes on 20 and 36 cells:
number and energy drift stay below 7 × 10⁻¹⁶, every sampled
occupation remains positive, and entropy increases by about 0.9%.
A separate direct continuum integral of the derivative of the
nonconserved fourth momentum moment agrees to 0.35% between its
two quadratures. The discrete moments differ from its finer value
by 26.0%, 8.9%, and 3.3% under successive grid/angle refinement.
Closure alone does not guarantee an accurate redistribution rate.

On the stored early highly occupied spectra the same *formal bare*
flux within the frozen mass ansatz has total absolute node change
of order 10¹⁵–10¹⁶ per particle per reference H at x=0.05, over
a billion times the dilute estimate. This is a warning against
inserting just Bose factors into a vacuum event rate and calling
the result a physical relaxation rate. The .06 result shifts by
a factor about 2.8 across the compared momentum grids despite
relative angle stability on each grid; sampling of the narrow
peak limits convergence. The result is not a spectral width or
an actual optical depth. A finite-occupation in-medium kernel
is required before treating temporal evolution of these snapshots
as a justified kinetic prediction.

~~~sh
python -W error early_pole_scalar_elastic_stencil_reference.py
python -W error early_pole_scalar_elastic_stencil_reference.py --cells 36 --angle 12
python -W error early_pole_scalar_elastic_weak_moment_reference.py
python -W error early_pole_scalar_elastic_snapshot_flux_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Manufactured and weak-moment witness](results/early-pole-scalar-elastic-stencil-2026-11-26.json);
[formal frozen-snapshot fluxes](results/early-pole-scalar-elastic-formal-flux-2026-11-26.json).
