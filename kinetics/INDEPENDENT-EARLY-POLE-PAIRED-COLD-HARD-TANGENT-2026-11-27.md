# A matched hard-state tangent resolves the spurious warm-shape ordering reversal

In the [shared-grid pole-tangent comparison](INDEPENDENT-EARLY-POLE-PAIRED-COLD-TANGENT-TIME-2026-11-27.md), subtracting almost identical *hard-state* cold rates still left the orientation of the warm-shape comparison unstable. Keep the exact same native final n216 snapshot, held scale and exposure 2e−6. In the cold packet term only, evaluate a symmetric hard-moment perturbation about the common baseline along the actually evolved small hard difference. Apply the resulting directional tangent at the original, unboosted difference to both the scalar and reciprocal hard deposits. Continue using the existing anchored pole and reciprocal pole tangent. The warm packets and the other joined channels are not replaced by this approximation.

| Joined macro steps and symmetric hard probe multiplier | Markers: cold feedback | Distributed q4: cold feedback | Relative decrease in feedback magnitude |
| --- | ---: | ---: | ---: |
| 6, 512 | −7.90743e−19 | −7.64060e−19 | 3.37% |
| 6, 4096 | −7.91278e−19 | −7.63783e−19 | 3.47% |
| 12, 4096 | −7.92054e−19 | −7.63038e−19 | 3.66% |

Each joined test applies the same 5e−15 elastic step acceptance. All terminal energy defects stay below 7e−16. Moving the directional probe by a factor of eight or halving the macro width does not reverse this ordering. The feedback itself remains only about −3.4e−12 of the much larger cold gain. The hard probe multiplier does not multiply any physical warm collision rates; it enlarges only a centered derivative measurement before returning to the evolved difference.

The optional diagnostic applies to a near-null difference and rejects a large hard displacement. It reproduces the order of magnitude of the earlier amplified-flux checks, but does not establish a continuum physical difference or a universal error bound. The default path does not use the probe and its prior decay-only fixture is byte-identical.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
a=results/early-pole-joined-fine-216-memory-2026-11-27.json.gz
p=early_pole_live_scalar_exterior_paired_difference_reference.py
# Repeat for (n,b) = (6,512), (6,4096), and (12,4096).
n=6; b=4096
python -W error "$p" "$a" --duration 2e-6 --steps "$n" --order 8 --joined \
  --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15 --hard-probe-boost "$b"
python -W error "$p" "$a" --duration 2e-6 --steps "$n" --order 8 --joined \
  --anchor-pole --linearize-cold-pole --elastic-energy-budget 5e-15 --hard-probe-boost "$b" --warm-order 4
```

Raw: `results/early-pole-paired-cold-hard-tangent-2026-11-27.json`.
