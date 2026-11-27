# Place one cold kernel sample at a moment-matched centroid

The [cold-mesh audit](INDEPENDENT-FINITE-MASS-QED-SHARED-POLAR-COLD-MESH-2026-11-27.md)
preserves the total archived scalar weight by integrating the
underlying occupation across every cold cell rather than using a
single sampled occupation. Keep those same exact cell weights and
examine only where to sample the bare kernel.

Compare the standard center in q³ with the weighted mean q and
the exponential of the weighted mean log q within each cell.
The weights are the archived f_s(q) q²/E_s(q), not smooth
equilibrium weights. Also try splitting cells at the two exact
charged aperture crossings and evaluating one linear centroid
in each resulting subcell. The scalar weight integral is
preserved in every test; the prior exact-weight center rule
reproduces independently.

On the 32-cell fixture, signed relative errors versus the fully
split reference are:

| Grid | q³ center | weighted mean q | weighted mean log q | 34-evaluation split mean q |
|---|---:|---:|---:|---:|
| Fine | −4.9366e−4 | −6.6126e−5 | −3.2085e−5 | −6.6731e−5 |
| Third | +3.0522e−4 | +3.0542e−5 | +7.0029e−6 | +4.5058e−5 |

Weighted logarithmic placement uses no additional kernel evaluations
and improves these two cold32 tests by factors about 15 and 44.
It is not the best rule for every resolution: on the third
64-cell check the weighted mean q gives a smaller error. Merely
splitting at the crossings is not sufficient for a universally
better one-node quadrature, since cancellation across neighboring
cells also changes.

This concerns an aperture-scale bare finite part and a frozen
test weight. It does not recommend changing a physical cold-mode
kinetic integrator without an independent conservation and
evolution check.

Replay from kinetics:
```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error finite_mass_qed_shared_polar_cold_centroid_reference.py
```
[Raw placement comparison](results/finite-mass-qed-shared-polar-cold-centroid-2026-11-27.json).
