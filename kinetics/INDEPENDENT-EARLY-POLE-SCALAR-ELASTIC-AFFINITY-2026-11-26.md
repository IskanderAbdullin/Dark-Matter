# Keep the sign of the high-occupancy Bose gain minus loss

The two gross products in a bosonic 2→2 collision both contain the
same quartic occupation term. In the stored early spectra the peak
occupation exceeds floating-point resolution for adding one.
A direct evaluation of the gross products can therefore erase the
actual signed cubic contribution entirely.

For positive occupations, write each gross product as the same
product of occupied factors multiplied by the exponential of the
sum of log(1+1/f) over the opposite pair. Compute that log with a
stable log-sum and subtract the exponential factors with expm1.
This also applies to fractional positive node weights with two
particles on each side. At a common Bose temperature the log is
affine in energy; exactly matched discrete energy then cancels
the occupation factor, up to floating-point roundoff.

A static witness at occupancies around 10¹⁶ gives zero from the
direct subtraction but a signed result of 2.2 × 10⁴⁹ from both the
stable affinity and the independently canceled cubic expression,
within 10⁻¹³ relative. A fractional-node thermal test preserves
the balance identity under the same high-occupation conditions.
This is an algebra and numerical reference only. A solver still
needs consistent weights for forward and reverse channels,
treatment of exactly empty cells, verified phase-space convergence,
and a physical in-medium scattering model.

~~~sh
python -W error early_pole_scalar_elastic_affinity_reference.py
~~~

[Static witnesses](results/early-pole-scalar-elastic-affinity-2026-11-26.json).
