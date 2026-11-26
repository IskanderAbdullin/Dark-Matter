# Metamorphic tracer check: one density versus two equal-mass distributions

Use exactly the positive finite prompt-plus-recoil shape from the
[production-era diagnostic](INDEPENDENT-FINITE-COOLING-PAIR-RECOIL-2026-11-26.md),
not a new coupled production model. The 86 eV scalar is a tracer
against the same dominant CDM. The mixed PSD is either entered
as one species or as two same-mass species with identical total
background density and the separate positive prompt and recoil
tables. The independent dense input audit compares the unscaled
joined \(f_0\) and the sum of the two split \(f_0\)'s. The largest
moment-weighted L1 difference over orders 0,1,2,4 is about
\(1.61\times10^{-7}\). The recoil count share is 0.18579, not an
independent dark-matter density.

For the linear collisionless hierarchy in a fixed metric, written
for \(\delta f\) rather than \(\delta f/f_0\), the source is linear
in the unperturbed distribution derivative. The energy-weighted
component contrasts must add for the same total spectrum if
evolution and integration are exact. Each approximation needs
its own test: a single fluid closure with an averaged sound
speed need not give the sum of two separately closed fluids.

For example, in the nonrelativistic limit let \(f\) be the hot
energy fraction, \(c_c^2,c_h^2\) the two adiabatic coefficients,
and \(D_c,D_h\) their individual density contrasts. The pressure
assigned to the total by one mean coefficient minus the sum
of the two pressure closures, divided by total density, is
\[
f(1-f)(c_h^2-c_c^2)(D_c-D_h).
\]
It need not vanish when most of the number is cold and most
of the squared momentum is in the warm recoil. This is a closure
mechanism, not proof of the origin of any given numerical deviation.

## Independent linear controls

First, the isotropic geodesic displacement characteristic applies
\(j_0(k\,d)\) to the same positive number measure. The actual joined
input agrees with the sum of the two number-weighted components
at the 1.13×10⁻⁸ level, and quadrature refinement is below 5×10⁻⁸.
This is free streaming with no perturbation history or gravitational
forcing, not a transfer function.

Second, on a deliberately simpler Newtonian matter-dominated clock,
start with zero perturbations and force all species with the same
constant-sign prescribed potential. At the final scale \(a_f\), the
density response normalized to the matching zero-dispersion response
is the mean of a \(j_0\) kernel with weight
\((1-\sqrt{a/a_f})\,da\). The displacement for comoving momentum
\(q T_{\gamma,0}\) is
\[
{2 q T_{\gamma,0}\over m H_m}
(a^{-1/2}-a_f^{-1/2}).
\]
This positive weighting also gives the required superposition.
For initial scale 0.001 and z=20 at k=10 h/Mpc, the independent
control gives 0.90147 prompt, 0.37520 recoil, 0.80369 joined.
A common physical speed cutoff removes less than \(5.5\times10^{-12}\)
of either number sample at the start. The other two starting times
also preserve additivity. This is **not** native CLASS evolution
or a calibrator for its physical perturbation amplitude.

## A pending native numerical check, not a reported transfer

The one-component native runs at 320/224 momentum nodes and
two smoothings give nearly the same default-closure ratio at k=10:
0.88609, 0.72334 and 0.52247 at z=0, 5 and 20. Agreement at
two momentum resolutions does **not** settle the grouping issue.
Compared to the separately computed prompt-only run, simple
superposition would require an inferred recoil-only contrast of
roughly +0.48, −0.19 and −0.90 at those probes. These values are
not results of a recoil-only calculation. Treat the native numbers
as a diagnostic, not as physical suppression of either scalar or
total matter, pending split and closure-disabled controls.
The exact pressure integrals show why this test is necessary:
at z=20, recoil carries about 92.9% of the background scalar
pressure despite carrying only 18.6% of the scalar energy.

## Reproduction

Regenerate the original prompt and mixed PSDs and their native INI
inputs using the documented exporters. Then generate the two-component
native input (example directories):

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_pair_split_tracer_reference.py \
  mix_native split_native prompt_native --eos EOS.dat \
  --pair results/finite-cooling-during-production-pair-shape-2026-11-26.json
~~~

The ballistic and prescribed-potential checks each take EOS.dat,
the unremapped mixed σ=0.14 PSD, the unremapped prompt σ=0.14 PSD,
and the split generator's finite_pair_only_ratio5_s14.dat, in that order.

After running the two-species INI, compare it to the one-species and
prompt reference outputs:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_pair_split_tracer_reference.py   mix_native split_native prompt_native   --joined-psd psd_mix/ratio5-scalar-pair-mix-sigma0.14.dat   --prompt-psd psd/ratio5-scalar-sigma0.14.dat
~~~
