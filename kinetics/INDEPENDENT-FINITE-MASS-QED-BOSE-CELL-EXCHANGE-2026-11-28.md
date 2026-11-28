# A conservative exact update for one reversible Bose cell pair

The [timelike pole coefficients](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-AFFINE-SOURCE-2026-11-28.md) separate fixed charged-bath factors A and B. If a later kinetic calculation evolves both outgoing scalar and incident photon populations, it needs a numerically safe exchange update. This test builds the exact map for **one reversible edge with held charged reservoirs**, not a fit to the actual collective occupation or a complete interaction network.

With capacities v_s,v_p>0, occupancies S,P and the event rate f=A P(1+S)−B S(1+P), one event adds a scalar to the weighted cell and removes a photon. Thus d(v_s S)/dt=f=−d(v_p P)/dt, and N=v_s S+v_p P is invariant. For the scalar fraction x=v_s S/N the evolution is

```text
dx/dt = a(1−x) − b x + c x(1−x)
a=A/v_p, b=B/v_s, c=(A−B) N/(v_s v_p).
```

The physical interior root is unique: the logarithm of the forward/backward ratio is strictly decreasing along the fixed-N segment. It selects a stable root r, with local slope −D where D²=(c+a−b)²+4ab. A Riccati solution gives the exact map. Taking e=exp(−Dt) and H=1+(c/D)(x_0−r), evaluate the new scalar fraction as a positive weighted average of x_0 and r,

```text
x(t) = [e x_0 + (1−e) H r] / [e + (1−e) H].
```

The implementation computes the photon fraction with the mirror expression and solves its own equilibrium root with a↔b and c→−c, rather than subtracting a population from one. The positive endpoint forms of H and stable quadratic roots avoid cancellation in the large-occupancy limit. Thermal equal-temperature factors are a fixed point. Away from it the factorized rate times its log ratio is nonnegative, and the exact map cannot cross the unique balance point.

Independent DOP853 integration in 60 non-equilibrium controls agrees to 1.39e−13 absolute. Pair count balance agrees to 2.23e−16 relative; two exact substeps agree with one to 3.56e−15 absolute. A shared thermal fixture agrees to 2.03e−14 relative. The rate/clock units here omit the physical coupling and remaining phase weights.

The 24 archived-shape controls use the eight fixed-angle A/B pairs and three illustrative capacity ratios. The accurate minority final photon occupation ranges from 2.26e−9 to 5.66e−7; it is computed in the model, not imposed as an archive photon population. Its log affinity closes to 7.11e−15. For the actual large scalar occupancies, subtracting two float64 scalar states can lose a transfer of about .30 in this test, despite the relative total count check passing. The model therefore retains a separate signed reaction extent from the smaller photon inventory. The test records scalar, photon and charged-reservoir energy changes E_s Δ, −ω Δ and (ω−E_s) Δ; their event ledger closes to 6.94e−18 absolute. A production integrator should carry a compensated or separate event ledger and evolve the charged bath rather than silently dropping such transfers.

Run `python -W error finite_mass_qed_bose_cell_exchange_reference.py` from `kinetics/` for the [JSON receipt](results/finite-mass-qed-bose-cell-exchange-2026-11-28.json). No angular or momentum integrated rate, matching to spacelike cuts, actual pole occupation, or cosmological solution is inferred from this isolated exchange.
