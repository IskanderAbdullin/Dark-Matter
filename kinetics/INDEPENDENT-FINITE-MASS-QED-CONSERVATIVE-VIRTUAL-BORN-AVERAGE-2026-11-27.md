# Signed Born average with number-and-energy-conserving internal occupations

Combine the [two-moment reconstruction](INDEPENDENT-FINITE-MASS-QED-CONSERVATIVE-VIRTUAL-RECONSTRUCTION-2026-11-27.md)
with the [archived signed Born weighting](INDEPENDENT-FINITE-MASS-QED-BORN-WEIGHTED-CHARGED-OUTER-2026-11-27.md).
Smooth only the internal loop occupation, with a requested 5%
transition width. The archived external forward, inverse, and signed
weights, cold scalar mode, and bare Born aperture remain unchanged.
The reconstruction preserves the number and energy of every internal
cell and stays within the Pauli bounds in the separate validation.

Calculate the two positive-frequency bare charged cuts at external
quadrature nodes, split at every external occupation edge and any
interior reconstruction joint. Independently refine the radial and
outer integration orders and reduce the fractional Born endpoint
exclusion from .002 to .0005. The new calculation reproduces the
original integrated external weights. The difference of forward and
inverse integrals equals the signed integral before normalization.
The largest radial component mesh discrepancy is below 1e−6 and
the largest outer change from 6 to 10 nodes is below 5.6e−5 in the
bare coefficient units.

| Grid | Signed step average | Signed simple linear average | Signed conservative average |
|---|---:|---:|---:|
| Fine | −1181.53906 | −1125.69888 | −1191.50769 |
| Third | −1559.35472 | −1374.65945 | −1584.95524 |

On the tighter support the conservative change increases the
magnitude by about .844% and 1.642%, respectively. This reverses
the direction of the ordinary nonconservative interpolation.
The forward and inverse answers respond similarly, and the quoted
relative signed changes shift by less than 3.4e−7 when tightening
the endpoint exclusion. This tests one requested smoothing family
and the selected archive snapshots; it does not bound an unknown
continuum occupation.

These are bare two-cut matrix coefficients with fixed external
kinematics. They are not corrections to the full production rate,
and do not replace physical widths, thermal cut matching, or
a phase-space update.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_conservative_virtual_born_average_reference.py
```
The replay includes independent inner and outer mesh checks and
can take multiple minutes. [Raw numerical receipt](results/finite-mass-qed-conservative-virtual-born-average-2026-11-27.json).
