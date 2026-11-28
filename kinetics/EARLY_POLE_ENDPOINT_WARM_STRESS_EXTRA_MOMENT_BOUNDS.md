# Extra even momentum moments can restore the sign for some heavy probes

The preceding constructive tests show that warm count and mean reference excess cannot determine the comparison sign at heavier formal masses. Here retain actual **additional per-cell even momentum moments** from each G3 and G4 nodal archive while keeping the same frozen cold convention.

For each warm cell put s=q² and take the center of its squared-momentum interval as s0. Write either kernel as a positive leading power of s times a binomial function of (s−s0)/(s0+M²). Fixing its first n+1 relevant even moments determines the corresponding degree-n expectation exactly. Bound the Taylor remainder analytically using the maximum absolute value of the in-cell scaled argument, and intersect it with a strictly positive bound from the leading moment and cell support alone. Check an independent direct nodal witness. These are ordinary floating-point enclosures, not directed-rounding certificates.

| Formal mass | Sufficient retained moments per warm cell for pressure sign | For pseudo-pressure sign |
| ---: | --- | --- |
| 1 | Not certified by this bound through degree 6 | Not certified by this bound through degree 6 |
| 3 | q², q⁴, q⁶, q⁸, q¹⁰ | Not certified through degree 6 |
| 10 | q², q⁴ | q⁴, q⁶, q⁸, q¹⁰ |
| 100 | q² | q⁴ |

At formal M=10 the two-moment pressure bound on G4/G3 minus one is [−0.48772%,−0.03521%]. At M=100 the one-leading-moment pressure bound is [−0.32313%,−0.22728%]; pseudo-pressure is [−0.85985%,−0.17136%]. Thus some higher moments are demonstrably useful, but “not certified” is no proof of ambiguity when those extra moments are actually retained. These masses remain formal probe parameters, not physical evolution or an observation.

From `kinetics/`:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=results/early-pole-fine-charge2-scalar4; d=2026-11-28
python -W error -m unittest test_early_pole_endpoint_warm_stress_extra_moment_bounds
python -W error early_pole_endpoint_warm_stress_extra_moment_bounds.py \
 --native "$p"-x05-state-$d.json.gz --gauss "$p"-gauss{3,4}warm-x05-state-$d.json.gz
```
