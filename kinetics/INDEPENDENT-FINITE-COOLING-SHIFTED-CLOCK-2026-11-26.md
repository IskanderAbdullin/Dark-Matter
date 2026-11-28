# Homogeneous clock sensitivity to the early converted population

The [retargeted spectrum](INDEPENDENT-FINITE-COOLING-SHIFTED-SPECTRUM-2026-11-26.md)
adds a positive early recoil spectrum at the same 1 MeV anchor.
Here remove the equal number of stable particles **without changing
the shape of their distribution**, then add that recoil spectrum.
The removed stable fractions are 0.05860% (instantaneous) and
0.09884% (finite). The prompt scalar count is unchanged, and
no later pair conversion is applied here.

For each recorded bin use its count and its mean p². Both the energy
and pressure integrands are increasing and concave in p². A chord
across a finite bin gives the lower estimate, evaluation at the
mean gives the upper, and the last interval's low edge supplies
a conservative lower estimate for that open tail. The two
background comparisons share the same stable shape. Subtract the
stable loss from the positive recoil before taking the square root,
rather than subtracting two unrelated total-density intervals.
These are bin-representation intervals within the stated proxy,
not uncertainties in the underlying production calculation.

| History | ΔH/H(rad) at 1 MeV | ΔH/H(rad) at 30 keV | Extra elapsed time from 1 MeV to 30 keV |
| --- | ---: | ---: | ---: |
| Instantaneous | [−1.1698, −1.1655]×10⁻¹⁰ | [−6.5366, −6.5361]×10⁻⁹ | about 6.36 μs |
| Finite, Tmax/TR=5 | [−2.3832, −2.3796]×10⁻¹⁰ | [−1.10661, −1.10657]×10⁻⁸ | about 10.79 μs |

The sign compares with the **unconverted** shifted prompt mixture,
not with a pure radiation cosmology. Removing rest mass decreases
density on this later interval even though the hotter recoil raises
pressure (ΔP/ρ(rad) about 3×10⁻¹¹ at these points). Numerically the
component proxies satisfy the collisionless continuity equation,
and an explicit moment-preserving positive discrete realization
lies inside the stated intervals at all four probes.

This does not calculate a self-consistent temperature response
during production, entropy injection, a nuclear abundance,
Delta N_eff, or a coupled nonlinear dark evolution. The tiny
clock shifts only measure this specific post-anchor comparison.

~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_cooling_shifted_collisionless_clock_reference.py \
  EOS.dat results/finite-cooling-shifted-{dist,spectrum-shape,count-audit}-2026-11-26.json
~~~
