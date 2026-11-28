# Charge refinement and a resolved cold-band test on evolved endpoint states

Two tests separate the evolution of the daughter spectrum from the quadrature of a particular conditional reaction. Both retain independent, opt-in archive-based diagnostics, not a coupled solution of the physical portal.

## Third-grid evolution with two charge subcells

Split each archived parent and fermion momentum cell into two before the empty-daughter startup, leaving the scalar grid fixed. Evolve the whole state from the initial clock through the four joins used by the earlier [third-grid endpoint replay](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-THIRD-CONTINUATION-2026-11-28.md). Restart each segment from the finer-time endpoint of the previous segment on the identical chart. The two fixed-start time meshes in each segment pass the stated pair-number, Pauli, and continuity tests.

| x | Third-grid two-subcell yield | Change from third-grid one subcell | Last fixed-start time change |
| --- | ---: | ---: | ---: |
| .01 | 3.47295534e−7 | −.034% | 1.31e−7 |
| .02 | 1.02151574e−5 | −.412% | 2.14e−6 |
| .04 | 1.71623040e−4 | −1.877% | 1.60e−4 |
| .05 | 3.65886183e−4 | −2.825% | 4.17e−5 |

At x=.05 the result is .493% above the original third archive, but 1.462% below the [fine two-subcell yield](INDEPENDENT-EARLY-POLE-ENDPOINT-FINE-CHARGE-REFINEMENT-2026-11-28.md). This is a larger fine-to-third spread than the .507% seen with one charge subcell. Refining charge within each grid does not by itself establish spatial convergence. In the third two-subcell card, the scalar pole mass is 8.68663043e−4, or .164% above the original third card, and the normalized scalar histogram total variation is .004463 from that original card. Neither is an uncertainty bound.

The [startup](results/early-pole-cell-endpoint-third-d2-startup-2026-11-28.json), [continuation](results/early-pole-cell-endpoint-third-d2-to-x05-2026-11-28.json) and [spectrum audit](results/early-pole-cell-endpoint-third-d2-spectrum-2026-11-28.json) preserve the numerical evidence. Reproduce using the prior startup driver with `third --div 2` and coarse step counts 2048 (empty start to .0075), 512 (to .01), then 1024 per segment through .05. Decompress the corresponding exact previous states from `results/early-pole-cell-endpoint-third-d2-*-state-2026-11-28.json.gz` for each `--restart`. The output state stores the fine-time trajectory for that segment.

## Cold aperture integrated on a single fine one-subcell state

Independently take the fine **one-subcell** x=.05 state, reconstruct both fermion occupations by conserving number and energy in each charge cell, and solve the transverse pole at fixed k=.0015. Split q in [1e−6,1.5e−5] at the scalar and charged support breaks. Keep the archived scalar shell convention and integrate all 27 panels using the smooth charged occupation callbacks on the exact radial support.

The [receipt](results/early-pole-fine-d1-endpoint-continuum-cold-band-2026-11-28.json) stores number, scalar energy and stable excess moments of the three nonnegative affine coefficients. The first two meshes differ by 1.573e−4 on the largest relative component and fail the predeclared 6e−5 gate; the next two agree to 2.173e−11. At this particular snapshot the integrated gross forward coefficient is 6.4584 times the [original fine result on the same band](INDEPENDENT-FINITE-MASS-QED-FIXED-POLE-COLD-SCALAR-BAND-2026-11-28.md), consistent with the earlier point probe but not inferred from it. The net scalar convention is `n_gamma D − B`, with `D=A−B`. The conditional number-neutral occupation is 3.34922719e−9; the scalar-energy-neutral occupation is larger by 4.47663e−8 **relative**, with the split changing by 2.20e−16 between the accepted meshes. This is neither a measured photon population nor the physical net production rate. The omitted common factor is `g² y² k² dk/(64π⁵)` for the same orientation and normalization used in the original band.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
gzip -dc results/early-pole-cell-endpoint-fine-x05-state-2026-11-28.json.gz > /tmp/fine-x05-state.json
for m in 0 1 2; do
  python -W error early_pole_endpoint_continuum_cold_band_pilot.py \
    /tmp/fine-x05-state.json "$m" "/tmp/fine-cold-$m.jsonl"
done
python -W error early_pole_endpoint_continuum_cold_band_receipt.py \
  results/fine-endpoint-continuum-cold-{1,2}.jsonl.gz
```

The two-subcell trajectories and the fine one-subcell band are distinct tests. The latter's resolved q integral does not make its evolved state spatially converged, and results for a full plasma and all photon and scalar momenta remain open.
