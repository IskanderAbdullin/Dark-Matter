# Preserve small entropy forces in manufactured high-Bose controls

The [unsplit averaged-entropy step](INDEPENDENT-TWO-MASS-DISCRETE-ENTROPY-GRADIENT-2026-11-27.md)
was tested at ordinary occupations. For very large Bose occupation
f=n/g, directly subtracting log(n+g) and log(n) can round the
entropy conjugate to zero even in extended precision. Evaluate
the two required positive line integrals directly instead:

    ᾱ_B = ∫₀¹ log1p(g/n(t)) dt,
    −∂ᾱ_B/∂n₁ = ∫₀¹ t g/[n(t)(n(t)+g)] dt,
    n(t)=n₀+t(n₁−n₀).

Split the quadrature uniformly in log count when the endpoints are
widely separated. This avoids subtracting two nearly equal
antiderivatives in the signal itself. The two line quantities and
the finite Bose entropy increment g_i s_B(f_i) differ from
independent 100-digit endpoint calculations by less than 2.23e−16
relative over the sampled two-endpoint fixtures. Their starting
occupations range from 10⁸ to 10³⁵ and count ratios from 10⁻⁸
to 10⁸. These large values are manufactured arithmetic controls.

There is a second potential cancellation when a tiny Bose exponent
is added to nearly identical parent and daughter Fermi exponents.
Keep the dominant and small pieces of each Fermi exponent separate.
Subtract the parent from each receiving fermion piece **before**
applying barycentric weights; then add the Bose contribution.
With normalized receiving weights that operation preserves the
physical event affinity but does not require recovering information
already lost in an absolute rounded input.

In three near-balanced manufactured events with Bose occupancy of
order 10²², the direct unsplit affinity rounds to zero. The split
evaluation resolves signs at imposed affinities ±10⁻²⁸ and 3e−29,
with signed flux magnitudes about 8.19e−7 and 2.46e−7 in the
fixture units. Independent direct 100-digit quantum products
confirm the results within 8.5e−14 relative.

These tests isolate the arithmetic. They do not demonstrate
convergence or stability of a full high-occupation nonlinear
solver, predict a condensate density, or supply physical
cosmological rates.

Replay from the kinetics directory:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error two_mass_saturated_bose_gradient_reference.py
python -W error two_mass_barycentric_small_affinity_reference.py
```
[Gradient receipt](results/two-mass-saturated-bose-gradient-2026-11-27.json)
and [affinity receipt](results/two-mass-barycentric-small-affinity-2026-11-27.json).
