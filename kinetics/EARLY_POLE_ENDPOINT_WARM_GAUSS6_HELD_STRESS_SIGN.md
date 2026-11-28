# Which G5/G6 stress signs follow from the actual cell moments?

The complete archived G5 and G6 snapshots at x=.05 retain the same cold and charge grid. Keep the cold cells at their original constant convention. In each warm cell constrain only its number and its mean excess energy at the common native S4 mass; do not assume the reconstructed exponential is the unique allowed profile.

At the actual measurement mass (about 8.73e−4), the independently certified relative G6/G5 pressure interval is [−0.000817662612, −0.000817503069]. Both pressure and pseudo-pressure also have a negative sign at every smaller sampled mass and at the sampled formal masses .01 and .1. At formal masses 1, 3, 10 and 100, the count-and-mean constraints allow **both** signs in both probes. Actual archived exponential witnesses are still negative; the sign is no longer determined by those cell moments alone.

The sign ambiguity is not an artifact of endpoint atoms. On each warm cell build bounded piecewise densities that are either centered or spread in normalized reference excess energy, with exact cell count and mean and cap .05 on each subinterval width. The same densities are used for pressure and pseudo-pressure. For example at formal mass 10 these give G6/G5 relative pressure −2.5303% or +2.1977%; the corresponding pseudo-pressure signs reverse together (−11.4905% or +11.8238%). Both signs occur with width cap .05 at all four heavier formal masses. Wider choices also admit both signs in both probes at masses 3, 10 and 100, but not necessarily at mass 1.

The analytically bounded extra-moment test starts with at least one additional moment. At formal mass 100, adding each warm cell’s actual nodal q² is sufficient for a negative pressure sign; q⁴ is sufficient for pseudo-pressure. At mass 10, the implemented bounds certify negative pressure after adding q² and q⁴, and negative pseudo-pressure after q⁴, q⁶, q⁸ and q¹⁰. These are sufficient under the tested hierarchy, not a necessity theorem. Formal observation masses do not map to physical times here.

```sh
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
d=2026-11-28; p=results/early-pole-fine-charge2-scalar4
python -W error early_pole_endpoint_warm_stress_conditional_hull.py --native "$p-x05-state-$d.json.gz" \
 --gauss "$p"-gauss{5,6}warm-x05-state-$d.json.gz \
 --receipts results/early-pole-frozen-stress-sweep-w124-gauss{5,6}-$d.json --levels 5 6
for n in early_pole_endpoint_warm_stress_extra_moment_bounds.py early_pole_endpoint_warm_bounded_shape_witness.py; do
 python -W error "$n" --native "$p-x05-state-$d.json.gz" --gauss "$p"-gauss{5,6}warm-x05-state-$d.json.gz --levels 5 6
done
```
