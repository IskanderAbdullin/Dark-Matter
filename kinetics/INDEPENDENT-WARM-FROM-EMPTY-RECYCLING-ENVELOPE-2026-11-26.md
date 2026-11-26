# Bound the full instantaneous 50 MeV vacuum charged-lepton benchmark

This independent follow-up tightens the kernel and start-to-late
energy/source envelopes of the
[project's recycling bound](https://github.com/IskanderAbdullin/Dark-Matter/blob/90c141722376dab64e4cd10422d545474817787c/kinetics/PHI-BOSE-PORTAL-RECYCLING-ENERGY-BOUND-2026-11-24.md).
It does not take the archived 10.442 MeV shapes as initial
conditions. Instead start with an empty dark sector at 50 MeV
after instantaneous reheating. Use only vacuum charged-lepton
pair input, its elastic heating, the project's specified
unpolarized dark pair network, and arbitrary internal momentum
redistribution. Bound the count even if Bose factors in the
allowed reverse/decay systems are large.

## Uniform energy-weighted kernel

The energy product inequality `E E′ ≥ s/4` upgrades the
CM rate to a lab rate per energy product. For the scalar
channel use the exact crossed vacuum matrix element, not a
sampling maximum. With `A=m₂²+s/2−m₁²−mφ²`,
`D²=4 p² k²`, `p²=s/4−m₁²`, `k²=s/4−mφ²`, and
`z=cos² θ`, the amplitude squared is
`8 y⁴ p² [s k⁴ z(1−z) + ((m₁+m₂)A−2m₁k²z)²] / (A²−D² z)²`.
Use positive monotone bounds in each angle/invariant cell.
The linear expression in the square increases with s and
decreases with z; use the largest squared corner, along
with the exact maximum of z(1−z) and a lower denominator
`(A²−D²)+D²(1−z)`. Above s=.8 GeV², use the analytic
[geometric denominator bound](INDEPENDENT-WARM-ISOLATED-RECYCLING-ENVELOPE-2026-11-26.md)
with β≤1 and the decreasing high-invariant coefficient.
The finite-s bounds decrease through 2.290, 2.141 and
2.071×10⁻¹² GeV⁻². The tail is at most 7.840×10⁻¹³.

For stable fermions, the opposite-charge angular integral
is at most `2(1+B/Δ)²`, from the same atanh inequality
with `B=(m₁+m₂)²−mφ²` and `Δ=m₁²−m₂²`. In the complete
low-invariant interval 4m₁²≤s≤.095 GeV², the concave
mixed-trace constant
`−s²/8+s(m₁²+m₂²+m₁m₂)−2m₁²m₂²` is positive at both
ends and hence throughout. The remaining angle-dependent
term is nonnegative, so the same-charge kernel is no
larger than the opposite-charge one there. The common
energy-weighted envelope peaks at s=5m₁². Above .095
GeV² the independent factor-two same-charge Cauchy bound
suffices because the envelope falls. Including two parents
per event and all incoming symmetry factors, the uniform
stable coefficient is ≤1.780×10⁻¹² GeV⁻² for any charge
abundances. Since the scalar and stable energies are disjoint,
the simultaneous coefficient is the larger of the scalar
and stable coefficients, not their sum. Sampled matrix
element folds are only independent diagnostics of these bounds.

## Positive source and cooling comparison

Enclose the nonzero-mass MB number and energy integrals
angle-independently in invariant-energy cells, using that
the phase factors and sub-pole propagator increase with
invariant mass below the split, while the Bessel functions
decrease. Use the full width floor and the exponentially
bounded tail at higher invariants. Include FD and final
Pauli by overestimating them with MB and zero blocking,
respectively. Positive comoving injection in each cooling
cell uses the rate at its warmer endpoint, which bounds the
whole cell by the positive invariant thermal kernel at fixed
invariant mass. Use the accumulated upper energy at the end of
each cell to enclose recycling there. The direct count includes
up to one decay per injected heavy, distinct from heavy pairs
regenerated internally. Gronwall retains the full elastic
heating allowance. The time tail uses a weaker analytic
massless injection and photon expansion ceiling when needed.

| Visible clock | Cooling cells | Direct upper contribution | Recycling upper contribution | Total Ωφh² upper |
| --- | ---: | ---: | ---: | ---: |
| Fixed radiation g*=10.75 | 250 | .00081385 | .04113035 | .04194420 |
| Fixed radiation g*=10.75 | 1500 | .00073387 | .03284822 | .03358209 |
| Tabulated central EOS, entropy conserved | 621 | .00059736 | .02280520 | .02340256 |
| Tabulated central EOS, entropy conserved | 1321 | .00057483 | .02096686 | .02154169 |

For the 50 MeV benchmark, both conditioned clocks yield an upper
scalar rest density below .060, the value corresponding to half
of .120. The more realistic central tabulated radiation clock
bounds the share of .120 by 18.0%. The finite-cell results are
upper comparisons, not best-fit estimates, and not observational
limits. If the total abundance is not .120, the ratio to .120
is not the actual scalar fraction.

The shorter late-anchor envelope and this from-empty upper
are different comparisons; do not stitch their outputs into
a trajectory. The result applies to these vacuum lepton
input and dark channels, not to late hadronic injection,
scalar number-changing vertices, hot-phase finite-width
plasma modifications or other reheating histories. The
numerical evaluation uses finite positive upper cells and
explicit analytic tails, with floating-point margins; it
does not claim machine-directed interval rounding. A relic
yield saturating the ceiling is not predicted.

~~~sh
export OPENBLAS_NUM_THREADS=1
D=results
python -W error warm_energy_weighted_recycling_kernels_reference.py > /tmp/kernels.json
python -W error warm_from_empty_recycling_envelope_reference.py /tmp/kernels.json --eos EOS.dat
~~~
