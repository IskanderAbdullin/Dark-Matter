# Doubling the fine charge resolution changes both the history and the cold aperture

Keep the fine scalar grid and the [cell-endpoint](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-COLLISION-2026-11-28.md) construction, but split each native parent and fermion bin into two intervals before the initial empty-daughter evolution. The original charge steps are retained at start. The horizon, solver gates and number/energy accounting match the earlier [one-interval trajectory](INDEPENDENT-EARLY-POLE-CELL-ENDPOINT-FINE-CONTINUATION-2026-11-28.md).

| x | Two-interval fine yield | Change from one interval |
| --- | ---: | ---: |
| .01 | 3.4754047e−7 | −.102% |
| .02 | 1.0219573e−5 | −.584% |
| .04 | 1.7314777e−4 | −2.404% |
| .05 | 3.7131297e−4 | −1.884% |

The x=.05 fixed-start time drift is 3.10e−5. The two-interval result is .303% above the old archived fine yield, but proximity is not a convergence bound. On the unchanged scalar grid its normalized histogram has total variation .01427 from that archive and the scalar pole mass is higher by .101%. The difference between one and two charge intervals grows through the history, despite their agreement much earlier at x=.0075. The next nested charge meshes must still be checked.

## A controlled change of the final cold aperture

At the two cold momenta q=3.5e−6 and 8.2e−6 with k=.0015, start with the corresponding original frozen archive. Change its scalar mass to the evolved value, then its photon pole to the moment-preserving continuum reconstruction, and finally both charged occupation profiles. This is a declared sequential replacement path, not an order-independent decomposition of a self-consistent plasma. The continuum occupation factors are integrated on the exact radial support.

At q=3.5e−6 the total raw forward factor relative to the respective original archive is 6.459 (fine one interval), 1.304 (fine two intervals), or 23.621 (third one interval). Along this path the mass alone gives 6.485, 1.314 and 23.621, respectively. The cosine-zero parent threshold changes strongly even when the scalar pole barely moves: for fine one interval it changes from 18.35145 to 16.76629 after both dispersion changes. The integrated occupation and actual photon population are not inferred by this point diagnostic. The cold raw response clearly requires stricter resolution than agreement of total scalar yield alone.

From `kinetics/`, decompress the needed `early-pole-cell-endpoint-*-state-2026-11-28.json.gz` states and run:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error early_pole_endpoint_cold_aperture_attribution.py FINE_X05.json FINE_D2_X05.json THIRD_X05.json
```
The [fine continuation](results/early-pole-cell-endpoint-fine-d2-to-x05-2026-11-28.json) retains all time segment reports; each segment replays via the earlier startup driver with `--div 2`. See the [spectral control](results/early-pole-cell-endpoint-fine-d2-spectrum-2026-11-28.json) and [conditional aperture receipt](results/early-pole-endpoint-cold-aperture-attribution-2026-11-28.json). The included third two-interval x=.0075 state is a restart, not its final history.
