# A rate-free 4-to-2 gap perturbation moves the narrow decay window

Differentiate the stationary Hartree gap with respect to changes
in the shell populations, holding the scale factor fixed. With
Eᵢ²=qᵢ²+m² and scalar bare mass m₀, the response to dnᵢ is

    dm = (λ/2) Σᵢ dnᵢ/Eᵢ / [2m + (λ/2)m Σᵢ nᵢ/Eᵢ³].

A cold threshold event changes total number by −2 and the
inverse-energy moment at fixed pole by −3/m. An initially
zero-momentum distribution therefore gives
d ln m/d ln n = 3(m²−m₀²)/[2(3m²−m₀²)], approaching one half
when the self-mass dominates. These statements do not assume a
rate or a population correlation closure.

On stored ratio-0.06 spectra, impose a chosen scalar number loss
by positive midpoint deposition at the current threshold, holding
both fermion distributions fixed. At 1% scalar number loss the
pole changes by about 0.5%, but the energy window for producing
a fermion with a partnered scalar at q=4.1 × 10⁻⁶ moves by
1.07 original half-widths at x=.05 and 1.93–1.94 at x=.10.
At x=.05 the pointwise average net stimulated factor moves from
0.0815 to 0.141 on the third grid, versus 0.0966 to 0.1013 on
the other grid. At x=.10 the relative changes agree better.
The geometry is more stable than the occupations within it;
no mesh-independent gain factor is inferred.

The independent tighter pole solve verifies ρ−3P equals
m₀²Σnᵢ/Eᵢ for the changing spectrum, including the populated hard
products. Midpoint count increments conserve total scalar
energy to better than 3 × 10⁻¹⁰ relative at the finer increment
count even for 5% imposed loss. The perturbation is a fixed-time
diagnostic, not a trajectory or a physical collision rate.

~~~sh
python -W error early_pole_chemical_gap_response_reference.py results/early-pole-spectrum-third-m06-2026-11-26.json.gz
~~~

[Numerical gap, trace and window checks](results/early-pole-chemical-gap-response-2026-11-26.json.gz).
