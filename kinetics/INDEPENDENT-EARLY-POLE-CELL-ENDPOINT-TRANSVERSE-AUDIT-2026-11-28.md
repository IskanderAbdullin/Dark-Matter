# A charge-cell quadrature can be adequate for evolution but not for the cold pole

The [cell endpoint construction](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-COLLISION-2026-11-28.md) reproduces each archived charge bin's number and first energy moment at its reference scale. Neither constraint fixes the different velocity-weighted moments in a plasma response. Recompute the collisionless transverse pole at k=.0015; compare direct endpoint summation, full archived constant-cell integration, and a continuous charge shape preserving the evolved moments.

Use the by-parts transverse response; for charge particle number measure dN its continuum kernel is proportional to `(2/E + v') H(βv) + β v v' H'(βv)`, with v=p/E, v'=m²/E³ and β=k/ω. This is finite at p=0, so cold endpoint nodes do not need a separate extrapolation. Independent polar orders and centered frequency derivative checks resolve the pole and residue. On the two original archive steps, the constant-cell reconstruction independently matches the previous response; the raw endpoint discrepancy tends down approximately fourfold per doubled charge subdivision.

| Original frozen archive | Raw node frequency change, 1 subdivision | Same at 8 | Corresponding change of cold forward coefficient, 1 subdivision | Same at 8 |
| --- | ---: | ---: | ---: | ---: |
| Fine | +2.71e−4 | +4.20e−6 | −18.14% | −.306% |
| Third | +1.58e−4 | +2.44e−6 | −12.65% | −.208% |

For that rate probe keep the original charge histogram and scalar mass fixed and change only the pole frequency/residue. It uses one cold q=3.5e−6, not an integrated correction. The high-energy aperture amplifies the optical discretization error even when the change in pole energy appears small.

For the [evolved fine x=.05 state](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-FINE-CONTINUATION-2026-11-28.md), a continuous Fermi-logit fit within each physical cell matches both charge number and first energy to 9e−14 in the normalized per-cell equations. With the preserved moments, its pole frequency differs by +2.90e−6 relative to the alternate number-preserving step fit. Raw endpoints differ from the step result by +2.75e−4. The reconstructed shape is a closure assumption, not an independently evolved continuum, and bounds on the actual pole or integrated rate do not follow from this comparison alone.

At q=3.5e−6 in the evolved fine snapshot, integrate that fitted occupancy inside the original radial panels using an opt-in callback, not by sampling it at panel midpoints. A constant-step callback reproduces the default branch to 2.23e−16. The continuous raw A, B and D agree between the two joint angular/radial orders to 7.75e−13 relative. Discretizing the same moment fit into 16 piecewise-constant subintervals leaves the inverse coefficient 3.67e−4 higher, even though its forward coefficient agrees much more closely. This is a reconstruction test at one point, not a rate integrated over q.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
gzip -dc results/early-pole-cell-endpoint-fine-x05-state-2026-11-28.json.gz > evolved-x05.json
python -W error early_pole_cell_endpoint_transverse_audit.py evolved-x05.json --rate
```
See the [receipt](results/early-pole-cell-endpoint-transverse-audit-2026-11-28.json). This is a separate frozen isotropic response, not a replacement of the default scalar or charge evolution.
