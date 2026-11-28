# At coincident charged poles the coherent bare residue is strongly suppressed

The [separated-root finite part](INDEPENDENT-FINITE-MASS-QED-CUT-LOCAL-SEQUENTIAL-SUBTRACTION-2026-11-27.md)
cannot be applied without checking the region where its pole
separation $h=2q\cdot K$ vanishes. On the spacelike Landau cut,
the exact locus is $\omega/k=q\cos\theta/E_s$ for positive
cosine. The spectral density is positive there, and both
charged roots coincide in occupied, Pauli-unblocked archive cells.
In the sampled cases their parent-energy locations agree to
2.23e−16.

Between external on-shell spinors the two projected numerators
are $N_1=2P\cdot\epsilon+\slashed K\slashed\epsilon$ and
$N_2=2b\cdot\epsilon-\slashed\epsilon\slashed K$.
Their difference is exactly $2q\cdot\epsilon$ times the scalar
spinor product. This follows from the gamma anticommutator
and $P-b=q-K$, without assuming a lightlike $K$.
When $D_2=-D_1$, the bare simple-pole numerator is the
**difference**, not the incoherent sum of the two squared
numerators. Across 48 occupied coalescence witnesses the
coherent difference-squared is only 6.98e−13 to 3.25e−12
of the incoherent sum. The direct spinor identity agrees
to better than 2.3e−11 relative to the individual numerators.

This result is a zero-width amplitude identity, not permission
to choose an imaginary sign or assign a rate at the merged pole.
For illustration, hold constant two real residue vectors with
this measured norm ratio $r$, normalize their summed diagonal
norm to one, and integrate over the full line. In units
$\pi/\epsilon$, using the same regulator sign in both
denominators gives integral one. Reversing the sign in just the
second gives
$(r+h^2/(4\epsilon^2))/(1+h^2/(4\epsilon^2))$.
Direct quadrature verifies both manufactured distributions.
This convention sensitivity when the separation is comparable
to the width is why neither independent diagonal subtraction
nor a zero-width coherent residue is, by itself, a physical
plasma rate. The retarded/advanced contour, sequential
overlap removal, statistical correlator, and compatible
dressing must be specified jointly.

Replay from kinetics:
~~~sh
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_cut_pole_coalescence_reference.py
OPENBLAS_NUM_THREADS=1 python -W error finite_mass_qed_cut_merged_width_convention_reference.py
~~~
[Occupied coherent roots](results/finite-mass-qed-cut-pole-coalescence-2026-11-27.json),
[width-sign algebra](results/finite-mass-qed-cut-merged-width-convention-2026-11-27.json).
