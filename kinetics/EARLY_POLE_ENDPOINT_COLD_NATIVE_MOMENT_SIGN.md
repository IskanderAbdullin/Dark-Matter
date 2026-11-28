# A second direct, sign-targeted split resolves the cold moment ordering

Start from the previously accepted eight-parent placement split. Use the
already measured panel counts to choose eight more midpoint splits
(four per grid). A balanced-child proxy proposes which splits to test.
The proxy does not determine the sign; evaluate every child directly
on both native quadrature meshes while retaining each chart's archived
kinetics, optical parameters and other certified rows. Join to
57 fine and 89 third panels on each mesh.

On the retained numerical q band, third-D2 is below fine-D2 for mean q,
RMS q, and mean speed at the fixed common reference mass. This holds for
forward number, inverse number and forward centered excess, under
arbitrary placement inside every final panel. For forward number the
third-minus-fine intervals are:

| Observable | Compatible interval |
| --- | ---: |
| Mean q | [−3.38460e−7, −1.40811e−8] |
| RMS q | [−3.45722e−7, −2.46608e−8] |
| Mean speed | [−3.87714e−4, −1.61297e−5] |

The two meshes agree on every sign. Reprojecting the child ratios onto
the exact original parent masses leaves all signs intact. The normalized
parent defects of the measured parent responses are below 1.81e−13 for
the three reported channel families on the reported mesh. This is a
conditional momentum ordering on the specified retained numerical band,
not a physical speed or a claim about the unmeasured support.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28
python -W error early_pole_endpoint_cold_native_moment_sign_targets.py --root results \
 --accepted "results/early-pole-cold-native-placement-eight-$d.json" \
 --archive "results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json"
python -W error early_pole_endpoint_cold_native_moment_sign_receipt.py --root results \
 --previous "results/early-pole-cold-native-placement-eight-$d.json" \
 --design "results/early-pole-cold-native-moment-sign-targets-$d.json" \
 --archive "results/finite-mass-qed-fixed-pole-cold-scalar-band-$d.json"
```
