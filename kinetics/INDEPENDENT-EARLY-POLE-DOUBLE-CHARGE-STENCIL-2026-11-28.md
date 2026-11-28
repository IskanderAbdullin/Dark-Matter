# A conserving alternative to a single daughter-cell centroid

The [axis test](INDEPENDENT-EARLY-POLE-ENDPOINT-AXIS-DESIGN-2026-11-28.md) isolates the daughter grid as the dominant source of susceptibility in the native endpoint decay clock. The native quadrature uses the aperture overlap of each daughter bin, but its parent deposition is tied to the energy of the entire daughter bin centroid, even for a partial overlap. Here keep the exact allowed energy interval for each scalar centroid and integrate within that interval, instead of evaluating the charged factors only at a daughter-bin representative.

Split the parent-energy aperture at the parent centroids and the shifted daughter centroids. Each event projects both charged momenta with nonnegative energy-interpolation weights on their respective neighbouring centroids. For each projected leg use a geometric average of its occupation factors and separately of its holes. The quotient of these averages has the exact common-temperature fermionic detailed-balance factor because the weighted centroid energies equal the event energies. Every event deposits one scalar and one daughter and removes one parent in the forward direction; the inverse is the same event with reversed net extent. The projection is separate from the scalar Hartree energy functional used by the unchanged snapshot.

On four combinations of native and number-remapped frozen states, relative event-energy errors are below 4.4e−16, thermal balance errors below 4.6e−15, and combined energy closure below 1.7e−15. Two quadrature orders agree on the gross forward, gross inverse and net clock to 1.45e−7 or better. Both tested charged cells have inward motion at zero and unit occupation. The restricted convex hull omits 1.50–2.27% of the thermal parent seed by this test's measure; that fraction is not a bound on the nonlinear archived net correction.

| Target mesh | Input archive | Native growth clock | Two-projection clock | Direct histogram clock on shared geometric support |
| --- | --- | ---: | ---: | ---: |
| Fine | Fine | 4.457997 | 3.359790 | 4.418265 |
| Fine | Third | 5.173667 | 4.452390 | 5.056491 |
| Third | Fine | 2.852987 | 3.146108 | 2.792807 |
| Third | Third | 2.189097 | 5.304473 | 3.670481 |

The new projection reduces the difference *between target grids for a fixed archived input*, but it does not uniformly reproduce the direct histogram integral. The latter samples exact step occupations on both charged grids and accumulates only the scalar source; it is not an energy-depositing three-species solver. Agreement between two destination meshes by itself therefore does not validate the modified collision source. Keep the alternative as a diagnostic until both the charged energy support and the representation of the sharp occupation steps are under control. No native integrator is replaced by this change.

Run `OPENBLAS_NUM_THREADS=1 python -W error early_pole_double_charge_stencil_pilot.py` in `kinetics/` and compare the [receipt](results/early-pole-double-charge-stencil-2026-11-28.json).
