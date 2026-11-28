# Preserve archived occupation steps in the two-charge projection

The [two-centroid stencil](INDEPENDENT-EARLY-POLE-DOUBLE-CHARGE-STENCIL-2026-11-28.md) conserves projected energies and respects thermal balance, but its wide interpolation across a sharp archived charge edge can change the instantaneous net rate. This test checks that approximation independently at the same two real frozen archive endpoints.

First subdivide the daughter bins uniformly, leaving parent and scalar grids fixed. An independent scalar integral on a shared energy interval splits at each native histogram edge, and can use exact step factors, parent-only geometric factors, daughter-only geometric factors, or both. At eight child cells per native daughter bin the daughter's interpolation error decreases, but the parent-only bias remains. For example, at that subdivision the fine direct, parent-only and two-interpolated clocks are 4.41827, 3.96692 and 3.94592; the third gives 3.67048, 2.80879 and 3.17362. Refining the daughter alone does not recover the direct reaction source.

Next retain each native charged bin as three subcells: two very thin ones against the faces and a central one. The archived step is unchanged; each subcell starts with the original constant occupation, keeping number and the Pauli bounds. Since all subcells within one old bin agree initially, the geometric transition is localized to a thin region across an actual archived edge. The event energies are still projected with two neighbouring centroids on both charged species. Unlike global uniform subdivision this uses exactly three times the charged bins per leg at every face-width refinement (and the same scalar mesh).

| Archive | Unpadded relative net error | Both padded, wing width 1/32768 of each cell | Largest gross coefficient error at that width |
| --- | ---: | ---: | ---: |
| Fine | −.239568 | +2.782e−6 | 7.368e−6 |
| Third | +.445171 | +2.074e−5 | 2.857e−5 |

Each comparison uses the independent direct piecewise-constant histogram integral over the **same** allowed energy interval. At the narrowest pad, parent-only or daughter-only refinements still have substantial net biases; the combined pad passes the stated 4e−5 coefficient check. The omitted thermal-parent-seed fraction drops to about 6.49e−12 on these frozen cases because the charged centroid hulls approach the physical edges. Exact event-energy, paired number, thermal detailed balance and selected Pauli boundary tests continue to hold. The cutoff fraction is a thermal seed diagnostic and not a bound on a general far-from-equilibrium tail.

This is an instantaneous numerical construction on fixed inputs. Evolving the thin subcells can separate their occupations and introduce a different stiffness and grid-convergence problem. No default driver is changed.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_nested_daughter_projection_pilot.py
python -W error early_pole_edge_padded_charged_projection_pilot.py
```
Compare the [nested daughter receipt](results/early-pole-nested-daughter-projection-2026-11-28.json) and the [edge-padding receipt](results/early-pole-edge-padded-charged-projection-2026-11-28.json).
