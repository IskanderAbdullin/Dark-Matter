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

## Native numerical grouping and matched closure-switch checks

The mixed scalar is an intentionally extremely subdominant tracer.
The default-fluid joined runs at 320 and 224 nodes agree with
each other to 1.78×10⁻⁵ across their common test range, but they do
not agree with a split into two positive same-mass tracers
having the same total \(f_0\). At k=10 h/Mpc:

| z | One mixed species, default closure | Two species, energy-weighted, default closure | Split prompt | Split recoil |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.886091 | 0.930709 | 0.978646 | 0.720632 |
| 5 | 0.723342 | 0.845855 | 0.932780 | 0.464919 |
| 20 | 0.522470 | 0.737279 | 0.847093 | 0.256034 |

These are amplitudes relative to CDM, not powers. The separated
prompt component agrees with its standalone run to better than
\(7\times10^{-10}\) in these probes. The 0.215 discrepancy at
z=20 far exceeds the input spectrum sum error and the change
under the two joined momentum quadratures. A negative recoil
contrast *inferred* earlier by subtracting the standalone prompt
from the joined-fluid result is therefore not a valid
recoil-only native calculation.

A separate matched six-mode test keeps the same joined input,
\(\ell_{\max}=36\) and exact k, and changes **only** the native
fluid approximation switch. At z=20:

| k [h/Mpc] | Fluid enabled | Fluid disabled |
| --- | ---: | ---: |
| 4 | 0.786468 | 0.881297 |
| 10 | 0.522605 | 0.759660 |

The two time histories agree to within 3.70×10⁻⁶ in scalar/CDM
amplitude on the sampled pre-switch intervals. At k=10 the raw
trajectory changes after the switch at conformal time 2.969 Mpc.
CDM contrasts at z=20 change by at most 1.3×10⁻⁷ relative.
This isolates a strong numerical fluid-closure dependence in
the reported tracer. It does **not** establish all-scale
convergence of either the switch-disabled finite hierarchy or
the separately closed split case. A converged grouped-vs-split
kinetic run is still required.

## Pressure mechanism check, not a calibration

At the native z=20 split contrasts held fixed, assigning the
joined adiabatic pressure coefficient gives 2.475 times the
pressure assigned by the two individual coefficients. This uses
the exact background integrals and is an algebraic evaluation
of the closures at one common state, not a claim that the
separately integrated pressure histories differ by that factor.

An independent controlled fluid ODE on the same simple EdS
constant-potential, zero-initial clock uses only
\(c_i^2(a)=5\langle q^2\rangle_i T_{\gamma,0}^2/(9m^2a^2)\).
It independently gives 0.56260 (one coefficient) versus 0.79332
(energy-weighted separate coefficients) at z=20 and k=10 with
initial \(a=0.001\). This is a mechanism demonstration, not a
CLASS fit. The independent **kinetic** EdS response remains
additive as above. Using the recorded same-sign native metric
as a late source in the kinetic control, also with zero initial
data at \(a=0.001\), yields 0.80525 joined and a sum defect less
than \(9\times10^{-9}\). Neither controlled zero-initial solution
is the native growing mode.

## Duplicate control and wider retained-fluid test

As a negative control, split a single mixed table into two
identical scalar species with exactly half the density each,
the same temperature and quadrature, and no change to its summed
background. Use the explicit RK evolver for both the single and
duplicated configs, with the default fluid approximation enabled.
At k=4 and k=10 h/Mpc at z=20 the two energy-weighted results agree
within 3.92×10⁻¹¹ absolute. The one-species RK run agrees with the
previous one-species NDF result within 5.5×10⁻¹⁰. Thus the large
difference above is not a generic two-species count or integrator
effect from this check; it needs the different grouped spectra.

The separate lmax=22, later-trigger=42, wider-k≤22 fluid run still
gives 0.522375 at z=20,k=10, within 9.5×10⁻⁵ of the original
0.522470 there. It produces a numerical half crossing near
10.677 h/Mpc at z=20 and 20.219 at z=5; do **not** interpret those
crossings physically. Changing the usual hierarchy and switch
time without removing the grouping closure did not eliminate
the tested issue.

The unmodified solver can generate the duplicate pair with its
explicit Runge–Kutta evolver (evolver=0), avoiding the stiff
evolver's dense Jacobian storage. This is a verified focused control,
not a general promise of lower wall time or convergence.

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

Generate the focused INIs from existing joined and split ones in
the native output directory, run the requested cases with the
native solver, retain its logs, and compute matched comparisons:

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --prepare native
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --cases joined_on_l36 joined_off_l36
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --time-probe
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --cases joined_on_l36 joined_on_l36_rk duplicate_on_l36_rk
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_focused_closure_reference.py native --wide
~~~

The mean-fluid control takes the same four arguments as the EdS
potential control. For the recorded-metric drive, add the focused
no-fluid native perturbation table for k1 after EOS.dat and before
the three PSD files.
