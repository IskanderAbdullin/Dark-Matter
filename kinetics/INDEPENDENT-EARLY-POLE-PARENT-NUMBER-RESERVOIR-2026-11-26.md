# Finite parent count matters more when scalar pole feedback is absent

The [paired daughter matrix](INDEPENDENT-EARLY-POLE-PAIRED-HARTREE-MATRIX-2026-11-26.md)
holds the parent FD occupation exactly.
Here keep its momentum **shape**
but allow its overall normalization $\alpha$
to reflect a finite evolving parent count:
$f_P=\alpha f_P^\mathrm{FD}$,
$\alpha=n_P/n_{P,\mathrm{FD}}$.
This multiplicative shape prescription is not
a Fermi distribution with an evolving chemical
potential, and some unsolved physical interaction
would have to maintain it. All retained phases
keep the parent occupation below the Pauli bound.

The shared reaction counts satisfy
$d(n_P+n_\phi)/d\ln a=0$ for no external
replenishment, and $n_\phi=n_b$.
Use the sum $n_P+n_\phi$ as a solver variable.
Add only the chosen comparison source
$\gamma(n_{P,\mathrm{FD}}-n_P)$ per $\ln a$,
with $\gamma=0,3,100$ in three arms.
The externally held bath temperature and
radiation clock are unchanged. No such relaxation
coupling is inferred for the real system.

| Scalar feedback and endpoint | Scalar yield with no replenishment | With $\gamma=3$ | With $\gamma=100$ |
|---|---:|---:|---:|
| Small positive quartic, $x=.05$ | .00043167 | .00043711 | .00044322 |
| Zero quartic, $x=.025$ | .00279899 | .00360423 | .00476241 |

With the positive quartic the finite no-replenishment
result is only 2.7% below the held-parent result
on the same coarse grid; the depleted parent is
still at .945 of the equilibrium count at the
endpoint. For the no-feedback arm the parent is
at .641 of the equilibrium count without replenishment,
and the reduced yield is about 42% below the held-parent
comparison. The different $x$ values are the same
chosen comparison endpoints as in the paired driver,
so they are not redshift-matched cosmological branches.
The computation verifies pair counts and the finite
reservoir invariant when $\gamma=0$.
It does not evolve the actual parent momentum distribution,
reheating energetics, thermal self-energies or the
physical mechanism of replenishment.

~~~sh
python -W error early_pole_parent_number_reservoir_reference.py
~~~

[Prescribed reservoir histories](results/early-pole-parent-number-reservoir-2026-11-26.json.gz).
