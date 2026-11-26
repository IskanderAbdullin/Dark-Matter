# Resolving production and collisional occupation of initially empty soft cells

This adds a numerically vacuum-aware option to the independent finite-support
two-moment evolution from an empty dark sector. The visible source, collision
matrix elements, dark quantum statistics in the collision operator, and
prescribed visible history are unchanged. The prior warm-cut result is retained.

## Numerical method

A collision step solves the same backward-Euler balances for both cell moments.
In a cold cell whose initial occupation is tiny, the residual relative to the
old density can be enormous even when the first implicit iterate is physically
close. At each Newton iteration scale each pair of number and normalized
kinetic-energy residuals by the greater of old and trial count in that cell.
Hold the weights fixed through that iteration's line search, use a positive
common-equilibrium trial on participating cells, and require both moment
residuals at the accepted solution below 8e-12 on that scale. Using the common
equilibrium as a trial adds no particles to the conserved input; both counts
and energy are checked independently after evolution.

A finite source quadrature can give exactly zero in very soft bins, because
direct production requires highly boosted parents that lie beyond the ordinary
thermal integration window. If that case is explicitly enabled, retain zero
in the physical initial moments rather than setting a floor. A cell that
does not appear in any collision geometric panel retains its input moments.
A geometrically participating cell, even if initially empty, is included in
the implicit system. This is a resolution-aware numerical treatment of the
stated quadrature, not proof that a continuum portal source is exactly zero.
The dimensionless cumulative below-cut statistic counts the interior of the
represented momentum interval only; an optional first finite-volume cell
starts at q=0. It is not a separate condensate degree of freedom.

## Checks and caveats

- Record the same mesh with production and redshift alone, prompt return
  only, and both prompt return and below-pole elasticity.
- Refine time step and compare mesh endpoints/cell widths separately from
  the global heat budget. Use absolute counts and report the low/high support.
- These results rely on common earlier assumptions: fixed visible EOS,
  vacuum dark final-state factors in the portal source, no visible inverse
  reactions, instantaneous and isotropic heavy decay, no high-momentum
  return from outside the numerical domain. A finite cell adjacent to zero
  does not resolve a zero-momentum singularity or establish a condensate.

## Current benchmark (visible 50 to 45 MeV)

With the lower endpoint at zero, the 12-step combined runs have final
scalar count fractions below q/T_R=0.01 of 1.359e-7, 1.363e-7, and
1.364e-7 for first positive boundaries 0.005, 0.002, and 0.0002.
The last two calculations begin the collision evolution with two and six
scalar shells left empty by the source quadrature. Both final count/energy
ledgers close to 6e-14 relative. At q/T_R<0.005 the two finer inner grids
give fractions about 1.697e-8. This agreement is for integrated counts, not
for a pointwise spectrum. With only prompt return on the zero-to-16 mesh,
the 24-step fraction below 0.01 is 5.35e-44, essentially the same as the
source/redshift comparison on the nearly identical [.005,16] support.
On the first-positive-boundary 0.005 grid, the 24-step combined value is 1.392e-7, about 2.45% above the 12-step value. The corresponding scalar collisional heat is 0.079011 of retained injected energy, compared with 0.079016 for prompt return alone.
The extreme source-only small number is merely the value under the stated
parent quadrature and not a converged precision prediction. Its contrast
with the nonlinear below-pole result identifies the channel populating these
soft bins under the model assumptions. The bulk scalar heat hardly differs
on the matched prompt-only and combined meshes, so a good bulk balance
alone does not test the cold population.

## Reproduce

From `kinetics`, using the same external central EOS as the previous checks:

```sh
export OPENBLAS_NUM_THREADS=1
BASE=(EOS.dat --physical-amplitude --visible-source --empty-pulse 1e-8)
SOFT=(--vacuum-scaled-residual --allow-empty-shells --include-zero-cell)
python -W error joint_resonant_frw_split_reference.py "${BASE[@]}" "${SOFT[@]}" --with-elastic --support-low .005 --empty-cells 18 --steps 12 24
python -W error joint_resonant_frw_split_reference.py "${BASE[@]}" "${SOFT[@]}" --with-elastic --support-low .002 --empty-cells 22 --steps 12
python -W error joint_resonant_frw_split_reference.py "${BASE[@]}" "${SOFT[@]}" --with-elastic --support-low .0002 --empty-cells 22 --steps 12
python -W error joint_resonant_frw_split_reference.py "${BASE[@]}" "${SOFT[@]}" --support-low .005 --empty-cells 18 --steps 12 24
python -W error joint_resonant_frw_split_reference.py "${BASE[@]}" --source-only --support-low .005 --steps 24
```
