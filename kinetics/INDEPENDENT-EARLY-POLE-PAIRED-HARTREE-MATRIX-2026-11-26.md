# Evolving the Pauli-blocked daughter reduces the held-parent early yield further

Extend the [shared scalar pole control](INDEPENDENT-EARLY-POLE-MULTIBAND-HARTREE-2026-11-26.md)
to evolve the fermion daughter in its own momentum
bins. Keep the parent FD and held, radiation clock
prescribed, narrow Lorentz-invariant masses, and
the same scalar Hartree-only closure. This is still
not a complete reheating or energy-closed model.

For each representative scalar momentum, the parent
on-shell energy interval is known analytically.
Intersect it with each daughter-fermion energy bin
shifted upward by the scalar energy. If the overlap
length is $D_{ij}$ and its parent-occupation
integral is $A_{ij}$, the signed collision
contribution before the common prefactor is
$$
 A_{ij}(1-f_{b,j})+
 f_{\phi,i}(A_{ij}-D_{ij} f_{b,j}).
$$
The same signed event-count matrix deposits both
species. The parent integral is analytic and uses
both finite endpoints. Rows with significant seed
match the independent held-daughter mode clock
to better than $10^{-10}$ relative in the test.
This caught an erroneous $x^2$ normalization
in the development control before production
results were retained.

The direct-occupation integrator preserves the
linear scalar-minus-daughter count invariant
to $4\times10^{-15}$ relative across the sampled
runs. A separate log-scalar coordinate integration
gives the same successful Hartree endpoint yields
but poorer pair-count conservation (up to $8.2\times10^{-6}$).
Its near-vacuum stiff control failed the accuracy
guard and was replaced by the direct integration.
All retained endpoints have fermion occupations
between zero and one, apart from allowed tiny
numerical tolerances. The matrix
conserves pair counts; bin-centroid energies are
not claimed to conserve total system energy.

| Reduced case | Evaluation $x=M/T$ | Scalar yield |
|---|---:|---:|
| Small quartic, previous empty fermion daughter | .05 | .00080224 |
| Small quartic, Pauli only, coarse | .05 | .00062006 |
| Small quartic, Pauli and inverse, coarse | .05 | .00044361 |
| Small quartic, Pauli and inverse, fine | .05 | .00043771 |
| Small quartic with a finite $b/T=.003$, coarse | .05 | .00049593 |
| Zero quartic with dynamical daughter, coarse | .025 | .00481969 |

The fine feedback endpoint is about 54.6% of the
empty-daughter comparison. On identical coarse
grids, allowing inverse reactions gives .715 of
the Pauli-only yield. That arm deliberately removes
the inverse process for comparison, not as a
thermodynamically complete model.
The coarse and fine massless-daughter feedback
yields differ by about 1.35%; the finite-daughter
and vacuum examples are single-grid diagnostics.
Nearly all of the coupled scalar yield in the
small-quartic case is still below $p/T=10^{-5}$
by momentum, not necessarily by stress.
No abundance was retuned after feedback.
An exact replenishing parent reservoir supplies
pair count and energy in this reduced operator;
neither its physical replenishment nor the
consistent thermal/chemical evolution is solved.

~~~sh
python -W error early_pole_paired_hartree_matrix_reference.py --direct
python -W error early_pole_paired_hartree_matrix_reference.py --direct --scalar 130 --fermion 190 --qlow 1e-9
python -W error early_pole_paired_hartree_matrix_reference.py --direct --no-inverse
python -W error early_pole_paired_hartree_crosscheck_reference.py "$RUN_DIRECTORY"
~~~

[Comparative ledger](results/early-pole-paired-hartree-crosscheck-2026-11-26.json)
and [fine direct spectra](results/early-pole-paired-hartree-fine-direct-2026-11-26.json.gz);
other runs have the same filename pattern with their
labels in the ledger.
