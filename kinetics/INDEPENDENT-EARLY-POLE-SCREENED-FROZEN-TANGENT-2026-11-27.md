# A thermal tangent does not describe the archived aperture perturbation

Linearize the finite screened elastic packet flux in the FD logits on the same prepared 97% parent states as the [short live witness comparison](INDEPENDENT-EARLY-POLE-NATIVE-WITNESS-JOINED-RELAXATION-2026-11-27.md). Hold all elastic geometry, screened weights and the energy gradients fixed at the reference state; do not include the derivative of the remap, active dielectric response or time drift of the baseline. Convert the logit derivative to a moment evolution using the FD susceptibility. Remove both charged-species number directions and total energy; for the same-flavor channel alone remove each flavor's energy separately.

As a numerical check, the projected matrix on the same held geometry at common affine-in-energy logits is symmetric and dissipative in the susceptibility metric. At the nonstationary archived logits it is not symmetric or uniformly contractive. A centered independently evaluated flux check matches the analytic derivative on the tested direction to better than 5e−8 relative, and the normalized conservation defects remain below 3e−13. The direction itself is the native-to-witness parent moment difference projected onto the appropriate constraint space; its logit derivative is scaled identically within each archived case for comparisons across sectors.

| Archive and held elastic operator | Fraction of initial susceptibility-weighted direction norm at exposure 1e−8, archived | Same fraction in affine control |
| --- | ---: | ---: |
| Fine, unlike only | .954 | .132 |
| Fine, same only | .667 | .129 |
| Fine, both | .659 | .0444 |
| Third, unlike only | .855 | .113 |
| Third, same only | .250 | .110 |
| Third, both | .240 | .0337 |

The same-flavor sector accounts for most of the initial decay of this chosen perturbation norm in the archived comparisons. In the both-sector fine case the norm rebounds from .578 at 3e−8 to .687 at 1e−7, consistent with the fact that there is no single thermal exponential describing this local model. A separate readout differentiates the scalar decay/inverse source with respect to the evolved perturbation. At 1e−8 this derivative is −.00228 on fine and +6.21e−5 on third; the third initial derivative was −6.19e−5. That is a sign change of a **directional derivative**, not of the net production source.

This diagnostic does not measure the true parent spectral width, prove a plasma instability or predict the nonlinear trajectory. The two archived states, meshes and prepared aperture directions differ; large sensitivity to those choices is itself part of the result. The equilibrium control is a manufactured comparison on the same finite graphs, not a calibrated thermal surrogate for the archived distributions.

From `kinetics/`, with one thread:
```sh
export OPENBLAS_NUM_THREADS=1
for g in fine third; do
 for s in both same unlike; do
  python -W error early_pole_screened_frozen_tangent_reference.py "$g" --sector "$s"
 done
done
```
Raw: `results/early-pole-screened-frozen-tangent-2026-11-27.json.gz`.
