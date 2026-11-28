# A redshift ledger below absolute-energy precision

Use the [kinetic-offset event projection](INDEPENDENT-TWO-MASS-KINETIC-OFFSET-TRACE-2026-11-27.md)
in a manufactured decay sequence with a prescribed scale factor.
Hold momentum-bin edges fixed in comoving coordinates. A unit initial
comoving parent count decays in four specified positive transactions
at a=1.2, 1.9, 3.6 and 6, with the final reading at a=8. Each event
uses a boost from the represented parent shell energy and is
positively deposited between neighboring daughter shell means.
The injection profile is not calculated from a collision rate.

Define the comoving kinetic energy K̃ and the comoving pressure P̃
by subtracting no rest masses in their numerical evaluation. For
each massive species let ΔΘ̃=mN−Θ̃ be the trace deficit. The direct
moments, independently formed from the same momentum distribution,
obey 3P̃=K̃+ΔΘ̃. Between events the exact redshift equation is
dK̃/dln a=−3P̃. Across one decay, the represented total kinetic energy
increases by the separately specified parent rest excess δ, while
the comoving rest energy decreases by the same amount. Thus

    K̃(a) − K̃(1) + ∫₁ᵃ 3P̃ dln a − δ N_decayed = 0.

The test evaluates both massive daughter trace deficits and the
parent contribution, plus pressure through an independent interval
quadrature. Both daughter counts match the removed parent count and
all kept occupations are nonnegative.

| Manufactured case | Largest kinetic-scale combined ledger defect | Largest per-event kinetic-scale defect |
|---|---:|---:|
| Ordinary rest excess and momenta | 1.38e−16 | 4.22e−18 |
| Near-rest, both absolute parent shell edges coincide | 5.12e−17 | 5.31e−18 |

At a=8 in the near-rest fixture, direct float64 subtraction of the
rest energy recovers zero kinetic energy. The split ledger still
resolves separate parent, scalar and fermion kinetic excesses of
1.27e−25, 3.03e−24 and 1.30e−24 in the receipt's units, together
with the pressure work of 5.70e−23. The precision figures refer
to this small kinetic budget, not its ratio to the much larger
rest mass.

This is a kinematic bookkeeping and regression reference. It does
not provide thermal masses, expansion history, source and inverse
reaction rates, or a physical relic abundance. The trace is not
imposed as a separate collision invariant.

Replay from the kinetics directory:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error two_mass_near_rest_redshift_ledger_reference.py
```
[Full receipt](results/two-mass-near-rest-redshift-ledger-2026-11-27.json).
