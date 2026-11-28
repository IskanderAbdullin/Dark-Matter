# Symmetric multi-edge Bose exchange with a retained reaction ledger

The [exact two-cell update](INDEPENDENT-FINITE-MASS-QED-BOSE-CELL-EXCHANGE-2026-11-28.md) preserves a pair's number and can retain its signed energy exchange when a scalar occupation is much larger than the photon occupation. Here one scalar cell shares reversible exchange edges with four photon cells. This is an operator-splitting check with the charged reservoirs and all edge coefficients held fixed, not a self-consistent collision solution or a way to infer a photon population.

Each edge has A P(1+S)−B S(1+P) in the convention of the prior test. Apply its exact two-cell map for half a step in one ordering, then for half a step in reverse ordering. Both submaps preserve nonnegative occupancies and total scalar-plus-photon count. Reserve an independent signed event extent per photon cell, computed from its inventory change; the net scalar increment is the sum of those event extents. Include the charged-reservoir energy with each event, because a scalar and a photon have different energies. A common Fermi/Bose bath temperature gives a fixed point of every edge, hence of the whole split map. The controllable test applies the appropriate entropy-affinity signs at each directed edge.

An independent DOP853 solve checks the coupled populations **and** the cumulative edge ledger. Halving the macro step lowers the maximum state-and-ledger error by factors 3.974, 3.993 and 3.998, consistent with second order on the nonstiff fixture. The shared thermal state is fixed to 2.23e−16 relative. The finite-step error itself is retained in the receipt and not claimed to vanish.

The two archived-shape stresses use the fixed-angle [timelike coefficient table](INDEPENDENT-FINITE-MASS-QED-TIMELIKE-POLE-AFFINE-SOURCE-2026-11-28.md), an actual cold scalar occupation, and illustrative photon inventories and capacities. They expose a bookkeeping failure:

| Archive | Scalar increase from signed photon extents | Repeated scalar float64 update | Resummed scalar float64 update |
| --- | ---: | ---: | ---: |
| fine | .1491266 | 0 | 0 |
| third | .1491089 | 12 | .25 |

On the third spectrum, repeated roundoff can build an apparent increase much larger than the real count. Reconstructing the stored scalar state from its original value and the signed event sum avoids that repeated drift, but the nearest scalar float64 value still cannot resolve the exact .149 increment. Thus keep the event ledger as a separate state quantity for number and heat rather than reporting a difference of scalar snapshots. With the event ledger, the energy closure residual is below 3.42e−20 in these normalized stress fixtures. No fitted timescale, equilibrium of the actual archive, full coupled plasma evolution, or matched scalar rate follows from the table.

Run `python -W error finite_mass_qed_bose_edge_network_reference.py` from `kinetics/` for the [machine receipt](results/finite-mass-qed-bose-edge-network-2026-11-28.json).
