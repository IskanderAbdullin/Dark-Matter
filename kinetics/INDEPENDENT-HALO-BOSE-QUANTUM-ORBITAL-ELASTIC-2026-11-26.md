# Final Bose factors enhance encounter counts more strongly in the central halo

Extend the [vacuum orbital diagnostic](INDEPENDENT-HALO-BOSE-ORBITAL-ELASTIC-CLOCK-2026-11-26.md)
by weighting each outgoing equal-mass contact
collision with both scalar final occupations.
The bath and orbit remain fixed. This is the
gross positive collision kernel, not the net
collision integral or a relaxation time.

For dimensionless pair CM speed $V$, relative
speed $g$, local Bose offset $y+\delta$, and
outgoing polar cosine $c$ relative to the CM
velocity, the two outgoing Bose exponents are
$a\pm bc$, where
$a=y+\delta+(v^2+u^2)/4$ and $b=Vg/2$.
Symmetry and the fixed sum $2a$ give the exact
polar average of the product $ (1+f_3)(1+f_4) $:
$$
 \frac{ b^{-1}\ln[(e^{a+b}-1)/(e^{a-b}-1)]-1}
      {1-e^{-2a}}.
$$
The continuous $b=0$ value is $(1-e^{-a})^{-2}$.
The formula is checked against independent numerical
polar integration. Both final factors are included;
no soft-angle or density-only multiplier is imposed.
Symmetry under exchange of the outgoing legs keeps
the angle-averaged tagged kick $g^2/2$ and energy
change $(u^2-v^2)/4$ per weighted encounter.

| Orbit in $R$, central occupation 7,000 | Gross stimulated encounters / vacuum |
|---|---:|
| Circular $.1R$ | 478 |
| Circular $R$ | 6.82 |
| Circular $10R$ | 1.048 |
| Circular $20R$ | 1.010 |
| $(.1,\,10)$ | 1.74 |
| $(.1,\,20)$ | 1.46 |
| $(5,\,20)$ | 1.07 |

The very central enhancement cannot be extrapolated
to an eccentric orbit that spends most of its time
outside the inner halo. For example the
$(.1,20)$ trajectory has 2.63 gross stimulated
encounters in 10 Gyr versus 1.80 without the Bose
factors, at the same illustrative cross section.
In the central toy trajectory the large count is
itself a warning against using an unperturbed orbit
as a dynamical solution. No claim about actual
halo evolution, condensation, or collisional
transport follows just from these counts.

~~~sh
python -W error halo_bose_quantum_orbital_elastic_reference.py
~~~

[Held orbital rate, kick and energy moments](results/halo-bose-quantum-orbital-elastic-2026-11-26.json).
