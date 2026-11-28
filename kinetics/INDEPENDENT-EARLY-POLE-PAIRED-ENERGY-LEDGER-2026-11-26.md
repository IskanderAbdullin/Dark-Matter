# Paired counts conserve exactly; bin-centroid energy does not

The [matched event counts](INDEPENDENT-EARLY-POLE-PAIRED-HARTREE-MATRIX-2026-11-26.md)
do not guarantee exact energy deposition into
representative bin momenta. This independent
diagnostic reconstructs the actual parent energy
absorbed from the held reservoir and compares it
against the energy assigned to the daughter bins.
For this test the fermion daughter is massless and
the scalar has the same repulsive self-consistent
pole used in the original paired feedback run.

For every overlap interval evaluate both
$\int f_P\,dE_P$ and the actual first moments
$\int E_P f_P\,dE_P$, $\int E_P\,dE_P$.
Compute the full **signed** gain minus reverse
energy, and keep a separate gross forward moment.
The projected daughter energy instead uses
$E_\phi+E_{b,j}$, with the representative
fermion bin energy. It need not equal the event
energy at points within that bin.
A separate continuity check uses the variational
Hartree energy, with its quartic subtraction, and
adds the explicit small vacuum-mass trace term
as $T\propto 1/a$. This trace is tiny on the
early interval but it is not set to zero.
Quadrature over the integrator's internal time
steps resolves the initial exponential growth.

This measures the microscopic *projection* error
before any correction of FRW bath evolution.
The total reservoir is held, so neither the true
finite thermalization rate nor a complete
stress-energy solution is inferred.

Both grids use a scalar cutoff of $10^{-9}T$ and end at $x=0.05$.
The signed continuum reservoir energy is 0.0236411 (coarse) and
0.0233796 (fine) in $T^4$ units. In the baseline centroid stencil,
the daughter energy assigned by collisions exceeds it by **0.5844%**
and **0.1854%**, respectively. The corresponding integral of the
*projected* source closes the daughter continuity equation to
2.8 × 10⁻⁷ relative on both grids, including the vacuum trace.
Thus a good continuity residual against a numerically defined source
does not detect the mismatch with the physical event moment.
The fine gross forward energy is 0.0353764; inverse decay removes
0.0119968 in the same units. The signed result is not the gross supply.

The [nodal control](INDEPENDENT-EARLY-POLE-PAIRED-NODAL-ENERGY-2026-11-26.md)
instead assigns a common energy to each discrete channel. Its energy
identity holds exactly in the discrete scheme, not against an exact
continuum FD integral.

~~~sh
python -W error early_pole_paired_energy_ledger_reference.py
python -W error early_pole_paired_energy_ledger_reference.py --scalar 130 --fermion 190
~~~

[Independent energy demand and projection checks](results/early-pole-paired-energy-ledger-2026-11-26.json).
