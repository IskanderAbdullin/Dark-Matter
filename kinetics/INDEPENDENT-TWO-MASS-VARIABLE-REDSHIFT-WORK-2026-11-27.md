# Redshift with prescribed varying masses needs an explicit trace-work ledger

Extend the [constant-mass near-rest redshift control](INDEPENDENT-TWO-MASS-NEAR-REST-REDSHIFT-LEDGER-2026-11-27.md).
Keep the comoving momentum grid fixed, but prescribe a different
smooth logarithmic mass slope for each daughter. Set the parent
mass to the sum of both daughters plus a separately stored positive
rest excess δ(a). Inject three known fractions of an initial parent
count; their rate is not computed.

For a single mode, with prime denoting d/dln a,

    E′ = −p²/E + (m′/m)(m²/E).
    K′ = −p²/E − (m′/m)(m−m²/E).

Thus the full one-particle energy receives work proportional to
the trace (m²/E), while the split kinetic equation uses its
deficit (m−m²/E). Summing over comoving populations, integrate
the two additional terms separately:

    ΔK̃ + ∫[3P̃ + Σ(m_i′/m_i)ΔΘ̃_i]dln a − Σ_events δ(a) ΔN = 0.
    Δρ̃_rest − ∫Σ m_i′ N_i dln a + Σ_events δ(a) ΔN = 0.

Both equations matter: the rest work from the prescribed mass
protocol is external to the particles. Assigning a new pole
mass without this work need not conserve their total energy.
The instantaneous decay transfers the rest excess to the kinetic
ledger with the opposite sign in the rest ledger; both daughter
count updates use positive shell weights.

| Manufactured case | Maximum kinetic-normalized defect | Maximum rest-normalized defect |
|---|---:|---:|
| Ordinary | 9.75e−17 | 4.72e−18 |
| Near rest | 7.40e−17 | 5.61e−18 |

In the near-rest example at a=8, the kinetic mass-deficit work is
−1.953e−24 in the fixture units, compared with positive pressure
work 5.904e−23. The separately integrated rest work is much larger
because it includes the rest energies of the prescribed particles.
The kinetic ledger resolves the small terms without subtracting
those nearly equal rest energies.

These are not self-consistent finite-temperature dispersion laws.
The source responsible for the mass change is deliberately external.
An energy-closed physical solver must derive the masses and account
for energy exchange with the medium, along with physical collision
and inverse rates.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_variable_redshift_work_reference.py
```
[Full receipt](results/two-mass-variable-redshift-work-2026-11-27.json).
